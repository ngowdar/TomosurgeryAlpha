using TomosurgeryAlpha;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Drawing;

namespace SpacerTests
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
        ///A test for LineSpacer
        ///</summary>
        [TestMethod()]
        public void LineSpacerTest()
        {
            float[,] f = new float[5,5]; // TODO: Initialize to an appropriate value
            RasterPath target = new RasterPath(f); // TODO: Initialize to an appropriate value
            int xstart = 10; // TODO: Initialize to an appropriate value
            int xend = 102; // TODO: Initialize to an appropriate value
            int[] expected = null; // TODO: Initialize to an appropriate value
            int[] actual;
            actual = target.LineSpacer(xstart, xend);
            Assert.AreEqual(expected, actual);
            Assert.Inconclusive("Verify the correctness of this test method.");
        }

        /// <summary>
        ///A test for ShotSpacer
        ///</summary>
        [TestMethod()]
        public void ShotSpacerTest()
        {
            float[,] f = null; // TODO: Initialize to an appropriate value
            RasterPath target = new RasterPath(f); // TODO: Initialize to an appropriate value
            int line = 0; // TODO: Initialize to an appropriate value
            PointF[] expected = null; // TODO: Initialize to an appropriate value
            PointF[] actual;
            actual = target.ShotSpacer(line);
            Assert.AreEqual(expected, actual);
            Assert.Inconclusive("Verify the correctness of this test method.");
        }
    }
}
