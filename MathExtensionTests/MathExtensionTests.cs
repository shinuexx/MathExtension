using Math = ShInUeXx.Numerics.Math;

namespace MathExtensionTests
{
    [TestClass]
    public class MathExtensionTests
    {
        [TestMethod]
        public void ValidateIfIsIntegerWorks()
        {
            var integer1 = 3d;
            var integer2 = 3f;
            var integer3 = 3m;

            Assert.IsTrue(Math.IsInteger(integer1));
            Assert.IsTrue(Math.IsInteger(integer2));
            Assert.IsTrue(Math.IsInteger(integer3));

            var non_integer1 = System.Math.PI;
            var non_integer2 = MathF.PI;
            var non_integer3 = 1.001m;

            Assert.IsFalse(Math.IsInteger(non_integer1));
            Assert.IsFalse(Math.IsInteger(non_integer2));
            Assert.IsFalse(Math.IsInteger(non_integer3));
        }
        [TestMethod]
        public void CanCallExternFunction()
        {
            Assert.AreEqual(0d, Math.Erf(0d));
        }
    }
}