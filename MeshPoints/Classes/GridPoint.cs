using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino.Geometry;

namespace MeshPoints.Classes
{
    public static class Globals
    {
        public const int GRID_RESOLUTION = 20;
        public const int PATCH_SIZE = 2;
        public const float VAL_RANGE = 2.86F;
    }

    public class GridPoint
    {
        public int PointId { get; set; }
        public double X { get; set; }
        public double Y { get; set; }
        public double Score { get; set; }
        public List<GridPoint> Neighbours { get; set; }

        public GridPoint(int _PointId, double _X, double _Y) {
            // What goes in a constructor? idk
            PointId = _PointId;
            X = _X;
            Y = _Y;
        }

        public static List<List<GridPoint>> GeneratePointGrid()
        {
            var pointGrid = new List<List<GridPoint>>();
            int pointId = 0;
            List<double> xCoordinates = new List<double>();
            List<double> yCoordinates = new List<double>();

            // Calculate x- and y-coordinate values
            int gridDimension = Globals.GRID_RESOLUTION + 2;
            for (int i = 0; i< gridDimension; i++)
            {
                double value = i / gridDimension * Globals.VAL_RANGE;
                xCoordinates.Add(value);
                yCoordinates.Add(value);
            }

            // Build empty point grid
            foreach (var y in yCoordinates)
            {
                List<GridPoint> pointGridRow = new List<GridPoint>();
                foreach (var x in xCoordinates)
                {
                    pointGridRow.Add(new GridPoint(pointId, x, y));
                    pointId += 1;
                }
                pointGrid.Add(pointGridRow);
            }

            // Set neighbours for each grid point (expect for the pads)
            for (int y = 1; y < gridDimension-1; y++)
            {
                for (int x = 1; x < gridDimension - 1; x++)
                {
                    List<GridPoint> neighbourList = new List<GridPoint>();
                    neighbourList.Add(pointGrid[y - 1][x - 1]);
                    neighbourList.Add(pointGrid[y - 1][x]);
                    neighbourList.Add(pointGrid[y - 1][x + 1]);
                    neighbourList.Add(pointGrid[y][x + 1]);
                    neighbourList.Add(pointGrid[y + 1][x + 1]);
                    neighbourList.Add(pointGrid[y + 1][x]);
                    neighbourList.Add(pointGrid[y + 1][x - 1]);
                    neighbourList.Add(pointGrid[y][x - 1]);
                    pointGrid[x][y].Neighbours = neighbourList;
                }
            }

            // Remove padding
            var cleanedPointGrid = new List<List<GridPoint>>();
            foreach (var row in pointGrid)
            {
                List<GridPoint> cleanedRow = new List<GridPoint>();
                foreach (var point in row)
                {
                    // Only add a grid point if it has neighbours
                    if (point.Neighbours.Any())
                    {
                        cleanedRow.Add(point);
                    }
                }
                // Only add a row if it contains points
                if (cleanedRow.Any())
                {
                    cleanedPointGrid.Add(cleanedRow);
                }
            }

            return cleanedPointGrid;
        }

        public static List<List<GridPoint>> GeneratePatches(List<List<GridPoint>> pointGrid)
        {
            List<List<GridPoint>> patches = new List<List<GridPoint>>();
            for (int row = 0; row < Globals.GRID_RESOLUTION; row += Globals.PATCH_SIZE)
            {
                for (int col = 0; col < Globals.GRID_RESOLUTION; col += Globals.PATCH_SIZE)
                {
                    var p1 = pointGrid[row][col];
                    var p2 = pointGrid[row][col +1];
                    var p3 = pointGrid[row+1][col];
                    var p4 = pointGrid[row+1][col+1];

                    patches.Add(new List<GridPoint>(new GridPoint[] { p1, p2, p3, p4 }));
                }
            }
            return patches;
        }
        
        public static List<Point3d> InterpolateNodesFromGridScore(
            List<List<GridPoint>> pointGrid, 
            int internalNodeCount)
        {
            List<Point3d> internalNodes = new List<Point3d>();

            // Flatten point grid to make querying simpler
            List<GridPoint> tmpFlatGrid = pointGrid.SelectMany(x => x).ToList();

            for (int inc = 0; inc<internalNodeCount; inc++)
            {
                List<GridPoint> flatGrid = tmpFlatGrid;
                // Find min value of grid and get the grid point with the minScore.
                // todo: this may introduce bugs.
                double minScore = flatGrid.Select(gp => gp.Score).Min();
                GridPoint minGridPoint = flatGrid.Where(gp => gp.Score == minScore).FirstOrDefault();

                // Build a neighbourhood
                List<GridPoint> neighbourhood = minGridPoint.Neighbours;
                neighbourhood.Add(minGridPoint);

                double totalScore = flatGrid.Sum(x => Math.Pow(x.Score, 4));
                List<double> weights = new List<double>();
                foreach (var point in neighbourhood)
                {
                    weights.Add((totalScore - Math.Pow(point.Score, 4)) / totalScore);
                }
                double totalWeight = weights.Sum();

                // Interpolation
                double intX = 0.0;
                double intY = 0.0;
                int i = 0;
                foreach (var point in neighbourhood)
                {
                    intX += weights[i] * point.X;
                    intY += weights[i] * point.X;
                }
                intX /= totalWeight;
                intY /= totalWeight;

                internalNodes.Add(new Point3d(intX, intY, 0));

                // Remove neighbourhood from the pointgrid
                List<int> exclusionPointsById = new List<int>();
                foreach (var point in neighbourhood)
                {
                    exclusionPointsById.Add(point.PointId);
                }
                tmpFlatGrid.Clear();
                foreach (var point in flatGrid)
                {
                    if (!exclusionPointsById.Contains(point.PointId))
                    {
                        tmpFlatGrid.Add(point);
                    }
                }
            }

            return internalNodes;
        }
    }
}
