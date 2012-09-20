using TomosurgeryAlpha;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;

namespace Calc3DSliceDose_Test
{
    
    
    /// <summary>
    ///This is a test class for RasterPathTest and is intended
    ///to contain all RasterPathTest Unit Tests
    ///</summary>
    [TestClass()]
    public class RasterPathTest
    {


        private TestContext testContextInstance;

        /// <summary>
        ///Gets or sets the test context which provides
        ///information about and functionality for the current test run.
        ///</summary>
        public TestContext TestContext
        {
            get
            {
                return testContextInstance;
            }
            set
            {
                testContextInstance = value;
            }
        }

        #region Additional test attributes
        // 
        //You can use the following additional attributes as you write your tests:
        //
        //Use ClassInitialize to run code before running the first test in the class
        //[ClassInitialize()]
        //public static void MyClassInitialize(TestContext testContext)
        //{
        //}
        //
        //Use ClassCleanup to run code after all tests in a class have run
        //[ClassCleanup()]
        //public static void MyClassCleanup()
        //{
        //}
        //
        //Use TestInitialize to run code before running each test
        //[TestInitialize()]
        //public void MyTestInitialize()
        //{
        //}
        //
        //Use TestCleanup to run code after each test has run
        //[TestCleanup()]
        //public void MyTestCleanup()
        //{
        //}
        //
        #endregion


        /// <summary>
        ///A test for Calculate_3D_SliceDose
        ///</summary>
        [TestMethod()]
        public void Calculate_3D_SliceDoseTest()
        {
            float[,] f = null; // TODO: Initialize to an appropriate value
            RasterPath target = new RasterPath(f); // TODO: Initialize to an appropriate value
            DoseKernel dk = null; // TODO: Initialize to an appropriate value
            int slicethickness = 0; // TODO: Initialize to an appropriate value
            string savepath = string.Empty; // TODO: Initialize to an appropriate value
            target.Calculate_3D_SliceDose(dk, slicethickness, savepath);
            Assert.Inconclusive("A method that does not return a value cannot be verified.");
        }

        /// <summary>
        ///A test for WriteToFile
        ///</summary>
        [TestMethod()]
        [DeploymentItem("TomosurgeryAlpha.exe")]
        public void WriteToFileTest()
        {
            PrivateObject param0 = null; // TODO: Initialize to an appropriate value
            RasterPath_Accessor target = new RasterPath_Accessor(param0); // TODO: Initialize to an appropriate value
            string savepath = string.Empty; // TODO: Initialize to an appropriate value
            float[] sd = null; // TODO: Initialize to an appropriate value
            target.WriteToFile(savepath, sd);
            Assert.Inconclusive("A method that does not return a value cannot be verified.");
        }

        /// <summary>
        ///A test for ReadDoseFromFile
        ///</summary>
        [TestMethod()]
        public void ReadDoseFromFileTest()
        {
            float[,] f = null; // TODO: Initialize to an appropriate value
            RasterPath target = new RasterPath(f); // TODO: Initialize to an appropriate value
            string loadpath = string.Empty; // TODO: Initialize to an appropriate value
            float[] expected = null; // TODO: Initialize to an appropriate value
            float[] actual;
            actual = target.ReadDoseFromFile(loadpath);
            Assert.AreEqual(expected, actual);
            Assert.Inconclusive("Verify the correctness of this test method.");
        }
    }
}
