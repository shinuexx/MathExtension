using System.Runtime.InteropServices;

namespace ShInUeXx.Numerics
{
    internal static class NativeFunctions
    {
        [DllImport("cpp-math")]
        internal static extern double erf(double x);

        [DllImport("cpp-math")]
        internal static extern double erfc(double x);
        [DllImport("cpp-math")]
        internal static extern double tgamma(double x);
        [DllImport("cpp-math")]
        internal static extern double lgamma(double x);

        [DllImport("cpp-math")]
        internal static extern float erff(float x);
        [DllImport("cpp-math")]
        internal static extern float erfcf(float x);
        [DllImport("cpp-math")]
        internal static extern float tgammaf(float x);
        [DllImport("cpp-math")]
        internal static extern float lgammaf(float x);

        [DllImport("cpp-math")]
        internal static extern double expm1(double x);
        [DllImport("cpp-math")]
        internal static extern double log1p(double x);

        [DllImport("cpp-math")]
        internal static extern float expm1f(float x);
        [DllImport("cpp-math")]
        internal static extern float log1pf(float x);

        [DllImport("cpp-math")]
        internal static extern double exp2(double x);
        [DllImport("cpp-math")]
        internal static extern float exp2f(float x);
    }
}
