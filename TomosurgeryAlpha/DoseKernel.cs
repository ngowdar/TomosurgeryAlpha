using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Collections;
using System.Runtime.Serialization.Formatters.Binary;
using System.Threading;
using System.Threading.Tasks;
using System.Reflection;

namespace TomosurgeryAlpha
{
    public class DoseKernel
    {
        public float[][] dose;
        public static int N;
        int[] sectorsizes;
        double resolution;
        double[] location;
        public float[,] midplane;
        public DoseKernelInfo DKI;

        public DoseKernel(string path)
        {
            string hpath = GetHeaderFileName(path);
            FileInfo fi = new FileInfo(path);

            DKI = new DoseKernelInfo();
            if (hpath != null && fi.Exists)
            {
                LoadHeader(hpath);
                LoadDose(path);
                FindMidplane();
                DKI.Size = dose.GetLength(0);
                DKI.Name = fi.Name;
                DKI.GetListboxInfo();
                SetDosemidplaneForOpt();
            }
        }

        public DoseKernel(int size)
        {
            DKI = new DoseKernelInfo();
            switch(size)
            {
                case 4:
                    Create4mmKernel();
                    DKI.Name = "4mmKernel";
                    break;
                case 8:
                    Create8mmKernel();
                    DKI.Name = "8mmKernel";
                    break;
                case 16:
                    Create16mmKernel();
                    DKI.Name = "16mmKernel";
                    break;
                default:
                    Create4mmKernel();
                    DKI.Name = "4mmKernel";
                    break;
            }
            FindMidplane();
            DKI.Size = dose.GetLength(0);
            DKI.GetListboxInfo();
            SetDosemidplaneForOpt();                
        }

        private void Create16mmKernel()
        {
            using (Stream stream = Assembly.GetExecutingAssembly().GetManifestResourceStream("16mm_head.bin"))
            using (StreamReader head = new StreamReader(stream))
            {
                LoadHeader(head);
            }
            using (Stream s = Assembly.GetExecutingAssembly().GetManifestResourceStream("16mm_dose.bin"))
            using (StreamReader t = new StreamReader(s))
            {
                LoadDose(t);
            }
        }

        private void Create8mmKernel()
        {
            string[] names = this.GetType().Assembly.GetManifestResourceNames();
            using (Stream stream = Assembly.GetExecutingAssembly().GetManifestResourceStream("TomosurgeryAlpha.Resources.8mm_head.bin"))
            using (StreamReader head = new StreamReader(stream))
            {
                LoadHeader(head);
            }
            using (Stream s = Assembly.GetExecutingAssembly().GetManifestResourceStream("TomosurgeryAlpha.Resources.8mm_dose.bin"))
            using (StreamReader t = new StreamReader(s))
            {
                LoadDose(t);
            }
        }

        private void Create4mmKernel()
        {
            string[] names = this.GetType().Assembly.GetManifestResourceNames();
            using (Stream stream = Assembly.GetExecutingAssembly().GetManifestResourceStream("TomosurgeryAlpha.Resources.4mm_head.bin"))
            using (StreamReader head = new StreamReader(stream))
            {
                LoadHeader(head);
            }
            using (Stream s = Assembly.GetExecutingAssembly().GetManifestResourceStream("TomosurgeryAlpha.Resources.4mm_dose.bin"))
            using (StreamReader t = new StreamReader(s))
            {
                LoadDose(t);
            }

        }

        public float[] ReturnDose1D()
        {
            float[] d = new float[N * N * N];
            for (int k = 0; k < N; k++)
                for (int j = 0; j < N; j++)
                    for (int i = 0; i < N; i++)
                        d[k * N * N + j * N + i] = dose[k][j * N + i];
            return d;
        }

        public float ReturnSpecificDoseValue(int i, int j, int k)
        {
            return dose[k][j * N + i];
        }

        private void FindMidplane()
        {
            float[] temp = dose[N / 2];
            midplane = new float[N, N];
            for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++)
                    midplane[i, j] = temp[j * N + i];
        }

        private string GetHeaderFileName(string path)
        {
            if (path.EndsWith("txt") == true)
            {
                path = Path.ChangeExtension(path, ".h");
                //path.Remove(count - 4,3);
                //String.Concat(path, "h");
                //path.Replace("txt", "h");
                //path.TrimEnd(("txt").ToCharArray());
                return path;
            }
            else
            {
                System.Windows.MessageBox.Show("You didn't choose a text file!");
                path = null;
            }
            return path;

        }

        private static ArrayList BinaryDeserialize(string path)
        {
            ArrayList output = null;
            using (FileStream str = File.OpenRead(path))
            {
                BinaryFormatter bf = new BinaryFormatter();
                output = (ArrayList)bf.Deserialize(str);
            }
            return output;
        }

        private void LoadDose(StreamReader t)
        {            
            dose = new float[N][];
            float[] temp;
            for (int j = 0; j < N; j++)
            {
                temp = new float[N * N];
                for (int i = 0; i < N * N; i++)
                {
                    temp[i] = (float)Convert.ToDecimal(t.ReadLine());
                    //dose[j][i] = (float)Convert.ToDecimal(t.ReadLine());
                }
                dose[j] = temp;
            }
        }

        private void LoadDose(string path)
        {
            StreamReader s = new StreamReader(path);
            LoadDose(s);
        }

        private void LoadHeader(StreamReader head)
        {            
            head.ReadLine();
            string sectorconfig = head.ReadLine();
            DKI.SectorConfig = Convert.ToInt32(sectorconfig);
            SetSectors(sectorconfig);
            head.ReadLine();
            resolution = Convert.ToDouble(head.ReadLine());
            DKI.Resolution = resolution;
            head.ReadLine();
            //resolution = Convert.ToDouble(res.Remove(0, 12));
            //head.ReadLine();
            location = new double[3];
            location[0] = Convert.ToDouble(head.ReadLine());
            location[1] = Convert.ToDouble(head.ReadLine());
            location[2] = Convert.ToDouble(head.ReadLine());
            DKI.ShotLocation = "" + location[0] + ", " + location[1] + ", " + location[2];
            N = Convert.ToInt16(head.ReadLine());
            DKI.Size = N;
        }

        private void LoadHeader(string hpath)
        {
            StreamReader head = new StreamReader(hpath);
            LoadHeader(head);
        }

        private void SetSectors(string sectorconfig)
        {
            sectorsizes = new int[8];
            for (int i = 0; i < 8; i++)
            {
                sectorsizes[i] = Convert.ToInt16(sectorconfig[i]);
            }
        }

        public float[] GetSlice(int z)
        {
            return dose[z];
        }

        public float[,] Get2DSlice(int z)
        {
            float[,] r = new float[N, N];
            float[] temp = dose[z];

            for (int j = 0; j < N; j++)
                for (int i = 0; i < N; i++)
                    r[i, j] = temp[i * j];
            return r;
        }

        public float[][,] GetDoseSlab(int startz, int endz)
        {
            float[][,] slab = new float[endz - startz][,];
            int s = startz;
            int e = endz;
            //Change start coords if requested start and end are beyond limits
            if (startz < 0)
                s = 0;            
            else if (startz >= N)
                s = -1;
            else
                s = startz;
            if (endz >= N)
                e = N - 1;
            else if (endz < 0)
                e = -1;
            else
                e = endz;

            if (s > 0 && e > 0)
            {
                Parallel.For(0, e - s, (i) =>
                {
                    slab[i] = Get2DSlice(i + s);
                });
            }
            return slab;

        }

        private void SetDosemidplaneForOpt()
        {
            float[,] dmp = new float[N, N];
            for (int j = 0; j < N; j++)
                for (int i = 0; i < N; i++)
                {
                    float z = 0;
                    for (int k = 0; k < N; k++)
                        z += dose[k][j * N + i];
                    dmp[i, j] = z;
                }
            RasterPath.dosemidplane = Matrix.Normalize(dmp);
        }

    }

    public struct DoseKernelInfo
    {
        public string Name { get; set; }
        public int Size { get; set; }
        public double Resolution { get; set; }
        public int SectorConfig { get; set; }
        public string ShotLocation { get; set; }
        public string Info { get; set; }

        public void GetListboxInfo()
        {
            Info = "DOSE: " + Name + ", " + Size + "x" + Size + ", Res: " + Resolution + "mm/px, Config: " + SectorConfig;            
        }
    }
}
