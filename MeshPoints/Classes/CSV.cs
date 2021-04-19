using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MeshPoints.Classes
{
    class CSV
    {
        /*Referer til: https://www.youtube.com/watch?v=vDpww7HsdnM */
        public static void addRecord(string variable, string filepath)
        {
            try
            {
                using (System.IO.StreamWriter file = new System.IO.StreamWriter(@filepath, true))
                {
                    file.WriteLine(variable);
                }
            }
            catch (Exception exeption)
            {
                throw new ApplicationException("Something went wrong.", exeption); 
            }
        }
    }
}
