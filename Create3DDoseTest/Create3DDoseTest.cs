using System;
using System.IO;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using TomosurgeryAlpha;


namespace Create3DDoseTest
{
    [TestClass]
    public class Create3DDoseTest
    {
        [TestMethod]
        public void TestMethod1()
        {
            //string s = "
            //DoseKernel dk = new DoseKernel();
            float[, ,] SliceSlab = GenerateTestSlab(5, 0);
            float[, ,] DoseSlab = GenerateTestSlab(3, 1);
            //float[][,] DoseSlab = dk.GetDoseSlab(startz, endz);
        }

        private float[, ,] GenerateTestSlab(int size, int value)
        {
            float[, ,] f = new float[size, size, size];
            for (int i = 0; i < size; i++)
                for (int j = 0; j < size; j++)
                    for (int k = 0; k < size; k++)
                        f[i, j, k] = value;
            return f;
        }

        private void LoadTestDK(string path)
        {
            //NEED TO FILL OUT
        }

        [TestMethod]
        public void TestMethod2()
        {

        }
    }
}
