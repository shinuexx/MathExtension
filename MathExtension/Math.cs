using System.Numerics;

[assembly: CLSCompliant(true)]
namespace ShInUeXx.Numerics
{
    /// <summary>
    /// Provides static extended static methods for mathematical functions.
    /// </summary>
    [CLSCompliant(true)]
    public static class Math
    {
        private static uint[] nums = { 1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1 };
        private static string[] rum = { "M", "CM", "D", "CD", "C", "XC", "L", "XL", "X", "IX", "V", "IV", "I" };
        private static Dictionary<char, uint> roman_dict = new()
        {
            {'I', 1},
            {'V', 5},
            {'X', 10},
            {'L', 50},
            {'C', 100},
            {'D', 500},
            {'M', 1000}
        };

        private const double rad_per_deg = System.Math.PI / 180;
        private const double deg_per_rad = 180 / System.Math.PI;

        private const float rad_per_deg_f = MathF.PI / 180f;
        private const float deg_per_rad_f = 180f / MathF.PI;

        private const double LN2 = 0.6931471805599453094172321214581765680755001343602552541206800094,
            A0 = 1.1975323115670912564578e0,
            A1 = 4.7072688112383978012285e1,
            A2 = 6.9706266534389598238465e2,
            A3 = 4.8548868893843886794648e3,
            A4 = 1.6235862515167575384252e4,
            A5 = 2.3782041382114385731252e4,
            A6 = 1.1819493347062294404278e4,
            A7 = 8.8709406962545514830200e2,

            B0 = 1.0000000000000000000e0,
            B1 = 4.2313330701600911252e1,
            B2 = 6.8718700749205790830e2,
            B3 = 5.3941960214247511077e3,
            B4 = 2.1213794301586595867e4,
            B5 = 3.9307895800092710610e4,
            B6 = 2.8729085735721942674e4,
            B7 = 5.2264952788528545610e3,

            C0 = 1.42343711074968357734e0,
            C1 = 4.63033784615654529590e0,
            C2 = 5.76949722146069140550e0,
            C3 = 3.64784832476320460504e0,
            C4 = 1.27045825245236838258e0,
            C5 = 2.41780725177450611770e-1,
            C6 = 2.27238449892691845833e-2,
            C7 = 7.74545014278341407640e-4,

            D0 = 1.4142135623730950488016887e0,
            D1 = 2.9036514445419946173133295e0,
            D2 = 2.3707661626024532365971225e0,
            D3 = 9.7547832001787427186894837e-1,
            D4 = 2.0945065210512749128288442e-1,
            D5 = 2.1494160384252876777097297e-2,
            D6 = 7.7441459065157709165577218e-4,
            D7 = 1.4859850019840355905497876e-9,

            E0 = 6.65790464350110377720e0,
            E1 = 5.46378491116411436990e0,
            E2 = 1.78482653991729133580e0,
            E3 = 2.96560571828504891230e-1,
            E4 = 2.65321895265761230930e-2,
            E5 = 1.24266094738807843860e-3,
            E6 = 2.71155556874348757815e-5,
            E7 = 2.01033439929228813265e-7,

            F0 = 1.414213562373095048801689e0,
            F1 = 8.482908416595164588112026e-1,
            F2 = 1.936480946950659106176712e-1,
            F3 = 2.103693768272068968719679e-2,
            F4 = 1.112800997078859844711555e-3,
            F5 = 2.611088405080593625138020e-5,
            F6 = 2.010321207683943062279931e-7,
            F7 = 2.891024605872965461538222e-15;

        /// <summary>
        /// Computes the <c>e</c> (Euler's number, 2.7182818) raised to the given power <paramref name="x"/>, minus 1.0. This function is more accurate than the expression <c><seealso cref="System.Math.Exp(double)">System.Math.Exp</seealso>(x) - 1d</c> if arg is close to zero.
        /// </summary>
        /// <param name="x">value of <see cref="double"/> type</param>
        /// <returns>If no errors occur <c>e^x - 1</c> is returned.</returns>
        public static double Expm1(double x) => InitializationClass.Initialize<double>() ?? NativeFunctions.expm1(x);
        /// <summary>
        /// Computes the natural (base <c>e</c>) logarithm of <c>1 + <paramref name="x"/></c>. This function is more precise than the expression <c><seealso cref="System.Math.Log(double)">System.Math.Log</seealso>(1 + <paramref name="x"/></c>) if <paramref name="x"/> is close to zero.
        /// </summary>
        /// <param name="x">value of <see cref="double"/> type</param>
        /// <returns>If no errors occur <c>ln(1 + <paramref name="x"/>)</c> is returned.</returns>
        public static double Log1p(double x) => InitializationClass.Initialize<double>() ?? NativeFunctions.log1p(x);
        /// <summary>
        /// Computes the <c>e</c> (Euler's number, 2.7182818) raised to the given power <paramref name="x"/>, minus 1.0. This function is more accurate than the expression <c><seealso cref="System.MathF.Exp(float)">System.MathF.Exp</seealso>(x) - 1d</c> if arg is close to zero.
        /// </summary>
        /// <param name="x">value of <see cref="float"/> type</param>
        /// <returns>If no errors occur <c>e^x - 1</c> is returned.</returns>
        public static float Expm1(float x) => InitializationClass.Initialize<float>() ?? NativeFunctions.expm1f(x);
        /// <summary>
        /// Computes the natural (base <c>e</c>) logarithm of <c>1 + <paramref name="x"/></c>. This function is more precise than the expression <c><seealso cref="System.MathF.Log(float)">System.MathF.Log</seealso>(1 + <paramref name="x"/></c>) if <paramref name="x"/> is close to zero.
        /// </summary>
        /// <param name="x">value of <see cref="float"/> type</param>
        /// <returns>If no errors occur <c>ln(1 + <paramref name="x"/>)</c> is returned.</returns>
        public static float Log1p(float x) => InitializationClass.Initialize<float>() ?? NativeFunctions.log1pf(x);

        /// <summary>
        /// Computes the <see href="https://en.wikipedia.org/wiki/Error_function">error function</see> of <paramref name="x"/>
        /// </summary>
        /// <param name="x">A value</param>
        /// <returns>Value of error function at <paramref name="x"/></returns>
        public static double Erf(double x) => InitializationClass.Initialize<double>() ?? NativeFunctions.erf(x);
        /// <summary>
        /// Computes the <see href="https://en.wikipedia.org/wiki/Complementary_error_function">complementary error</see> function of <paramref name="x"/>, that is <c>1 - <see cref="Erf(double)">Erf</see>(<paramref name="x"/>)</c>, but without loss of precision for large <paramref name="x"/>
        /// </summary>
        /// <param name="x">A value</param>
        /// <returns>Value of complementary error function at <paramref name="x"/></returns>
        public static double Erfc(double x) => InitializationClass.Initialize<double>() ?? NativeFunctions.erfc(x);
        /// <summary>
        /// Computes the <see href="https://en.wikipedia.org/wiki/Gamma_function">gamma function</see> of <paramref name="x"/>
        /// </summary>
        /// <param name="x">A value</param>
        /// <returns>The value of gamma function at <paramref name="x"/></returns>
        public static double TGamma(double x) => InitializationClass.Initialize<double>() ?? NativeFunctions.tgamma(x);
        /// <summary>
        /// Computes the natural logarithm of the absolute value of the <see href="https://en.wikipedia.org/wiki/Gamma_function">gamma function</see> of <paramref name="x"/>
        /// </summary>
        /// <param name="x">A value</param>
        /// <returns>The value of the logarithm of the gamma function at <paramref name="x"/></returns>
        public static double LGamma(double x) => InitializationClass.Initialize<double>() ?? NativeFunctions.lgamma(x);

        /// <summary>
        /// Computes the <seealso href="https://en.wikipedia.org/wiki/Error_function">error function</seealso> of <paramref name="x"/>
        /// </summary>
        /// <param name="x">value</param>
        /// <returns>Value of Error Function at <paramref name="x"/></returns>
        public static float Erf(float x) => InitializationClass.Initialize<float>() ?? NativeFunctions.erff(x);
        /// <summary>
        /// Computes the <seealso href="https://en.wikipedia.org/wiki/Complementary_error_function">complementary error</seealso> function of <paramref name="x"/>, that is <c>1 - <see cref="Erf(float)">Erf</see>(<paramref name="x"/>)</c>, but without loss of precision for large <paramref name="x"/>
        /// </summary>
        /// <param name="x">A value</param>
        /// <returns>Value of complementary error function at <paramref name="x"/></returns>
        public static float Erfc(float x) => InitializationClass.Initialize<float>() ?? NativeFunctions.erfcf(x);
        /// <summary>
        /// Computes the <see href="https://en.wikipedia.org/wiki/Gamma_function">gamma function</see> of <paramref name="x"/>
        /// </summary>
        /// <param name="x">A value</param>
        /// <returns>The value of gamma function at <paramref name="x"/></returns>
        public static float TGamma(float x) => InitializationClass.Initialize<float>() ?? NativeFunctions.tgammaf(x);
        /// <summary>
        /// Computes the natural logarithm of the absolute value of the <see href="https://en.wikipedia.org/wiki/Gamma_function">gamma function</see> of <paramref name="x"/>
        /// </summary>
        /// <param name="x">A value</param>
        /// <returns>The value of the logarithm of the gamma function at <paramref name="x"/></returns>
        public static float LGamma(float x) => InitializationClass.Initialize<float>() ?? NativeFunctions.lgammaf(x);

        /// <summary>
        /// Computes 2 raised to the given power <paramref name="x"/>
        /// </summary>
        /// <param name="x">value of <see cref="double"/> type</param>
        /// <returns>If no errors occur, the base-2 exponential of <paramref name="x"/> (<c>2^<paramref name="x"/>)</c> is returned.</returns>
        public static double Exp2(double x) => InitializationClass.Initialize<double>() ?? NativeFunctions.exp2(x);
        /// <summary>
        /// Computes 2 raised to the given power <paramref name="x"/>
        /// </summary>
        /// <param name="x">value of <see cref="float"/> type</param>
        /// <returns>If no errors occur, the base-2 exponential of <paramref name="x"/> (<c>2^<paramref name="x"/>)</c> is returned.</returns>
        public static float Exp2(float x) => InitializationClass.Initialize<float>() ?? NativeFunctions.exp2f(x);

        /// <summary>
        /// Converts the number in degrees to the radians equivalent
        /// </summary>
        /// <param name="degrees">Angular value in degrees</param>
        /// <returns>The radian equivalent of <paramref name="degrees"/></returns>
        public static double Deg2Rad(double degrees) => degrees * rad_per_deg;
        /// <summary>
        /// Converts the radian number to the equivalent number in degrees
        /// </summary>
        /// <param name="radians">A radian value</param>
        /// <returns>The equivalent of <paramref name="radians"/> in degrees</returns>
        public static double Rad2Deg(double radians) => radians * deg_per_rad;
        /// <summary>
        /// Converts the number in degrees to the radians equivalent
        /// </summary>
        /// <param name="degrees">Angular value in degrees</param>
        /// <returns>The radian equivalent of <paramref name="degrees"/></returns>
        public static float Deg2Rad(float degrees) => degrees * rad_per_deg_f;
        /// <summary>
        /// Converts the radian number to the equivalent number in degrees
        /// </summary>
        /// <param name="radians">A radian value</param>
        /// <returns>The equivalent of <paramref name="radians"/> in degrees</returns>
        public static float Rad2Deg(float radians) => radians * deg_per_rad_f;

        /// <summary>
        /// Normalized <see cref="Sinc(double)">Sinc</see> function. 
        /// </summary>
        /// <param name="x">A value</param>
        /// <returns>1 for <c><paramref name="x"/> == 0</c>; otherwise <c><see cref="System.Math.Sin">System.Math.Sin</see>(<paramref name="x"/>) / <paramref name="x"/></c></returns>
        public static double Sinc(double x) => x == 0 ? 1d : System.Math.Sin(x) / x;
        /// <summary>
        /// Normalized <see cref="Sinc(double)">Sinc</see> function. 
        /// </summary>
        /// <param name="x">A value</param>
        /// <returns>1 for <c><paramref name="x"/> == 0</c>; otherwise <c><see cref="System.MathF.Sin">System.MathF.Sin</see>(<paramref name="x"/>) / <paramref name="x"/></c></returns>
        public static float Sinc(float x) => x == 0 ? 1f : System.MathF.Sin(x) / x;
        /// <summary>
        /// Normalized <see cref="Sinc(double)">Sinc</see> function. 
        /// </summary>
        /// <param name="x">A value</param>
        /// <returns>1 for <c><paramref name="x"/> == 0</c>; otherwise <c><see cref="Complex.Sin">Complex.Sin</see>(<paramref name="x"/>) / <paramref name="x"/></c></returns>
        public static Complex Sinc(Complex x) => x == 0d ? 1d : Complex.Sin(x) / x;
        /// <summary>
        /// Calculate inverse of error function.
        /// </summary>
        /// <param name="x">A value</param>
        /// <returns>The value that satisfy <paramref name="x"/> equals <see cref="Erf(double)">Erf</see>(<see cref="ErfInv(double)">ErfInv</see>(<paramref name="x"/>))</returns>
        public static double ErfInv(double x)
        {
            if (x < -1 || x > 1) return double.NaN;
            else if (x == 1) return double.PositiveInfinity;
            else if (x == -1) return double.NegativeInfinity;

            var abs_x = System.Math.Abs(x);
            double r, num, den;
            if (abs_x <= 0.85)
            {
                r = 0.180625 - 0.25 * x * x;
                num = (((((((A7 * r + A6) * r + A5) * r + A4) * r + A3) * r + A2) * r + A1) * r + A0);
                den = ((((((B7 * r + B6) * r + B5) * r + B4) * r + B3) * r + B2) * r + B1) * r + B0;
                return x * num / den;
            }

            r = System.Math.Sqrt(LN2 - System.Math.Log(1 - abs_x));
            if (r <= 5)
            {
                r -= 1.6;
                num = (((((((C7 * r + C6) * r + C5) * r + C4) * r + C3) * r + C2) * r + C1) * r + C0);
                den = (((((((D7 * r + D6) * r + D5) * r + D4) * r + D3) * r + D2) * r + D1) * r + D0);
            }
            else
            {
                r -= 5;
                num = (((((((E7 * r + E6) * r + E5) * r + E4) * r + E3) * r + E2) * r + E1) * r + E0);
                den = (((((((F7 * r + F6) * r + F5) * r + F4) * r + F3) * r + F2) * r + F1) * r + F0);
            }

            return x < 0 ? -num / den : num / den;
        }
        /// <summary>
        /// Indicates wheater the value is integer.
        /// </summary>
        /// <param name="x">A value</param>
        /// <returns><see langword="true"/> if <paramref name="x"/> is integer; otherwise <see langword="false"/>.</returns>
        public static bool IsInteger(double x) => double.IsFinite(x) && ((x % 1d) == 0d);
        /// <summary>
        /// Indicates wheater the value is integer.
        /// </summary>
        /// <param name="x">A value</param>
        /// <returns><see langword="true"/> if <paramref name="x"/> is integer; otherwise <see langword="false"/>.</returns>
        public static bool IsInteger(float x) => float.IsFinite(x) && ((x % 1f) == 0f);
        /// <summary>
        /// Indicates wheater the value is integer.
        /// </summary>
        /// <param name="x">A value</param>
        /// <returns><see langword="true"/> if <paramref name="x"/> is integer; otherwise <see langword="false"/>.</returns>
        public static bool IsInteger(decimal x) => ((x % 1m) == 0m);
        /// <summary>
        /// Computes an approximation of the principal <paramref name="n"/>-th root of <see cref="BigInteger"/> as the largest integer less than or equal to <c>R</c> for which <c>Pow(R, <paramref name="n"/>)</c> == <c><paramref name="base"/></c>.
        /// </summary>
        /// <param name="base">A positive integer</param>
        /// <param name="n">A non-negative integer</param>
        /// <returns>A value that approximate <c>R</c> which satisfy <c>Pow(R, <paramref name="n"/>)</c> == <c><paramref name="base"/></c></returns>
        /// <exception cref="ArgumentOutOfRangeException">If <paramref name="n"/> is less than 1 or <paramref name="base"/> is less than 0.</exception>
        public static BigInteger Root(BigInteger @base, int n)
        {
            if (n < 1) throw new ArgumentOutOfRangeException(nameof(n));
            if (@base < 0) throw new ArgumentOutOfRangeException(nameof(@base));

            int n1 = n - 1;
            BigInteger n2 = n;
            BigInteger n3 = n1;
            BigInteger c = 1;
            BigInteger d = (n3 + @base) / n2;
            BigInteger e = ((n3 * d) + (@base / BigInteger.Pow(d, n1))) / n1;
            while (c != d && c != e)
            {
                c = d;
                d = e;
                e = ((n3 * e) + @base / BigInteger.Pow(e, n1)) / n2;
            }
            return d < e ? d : e;
        }
        /// <summary>
        /// Returns the angle whose hyperbolic cosine is the specified <see cref="Complex"/> number.
        /// </summary>
        /// <param name="value">A <see cref="Complex"/> number representing a hyporbolic cosine.</param>
        /// <returns>An angle, mesured in radians, which is the arc hyperbolic cosine of <paramref name="value"/></returns>
        public static Complex Acosh(Complex value)
        {
            var a = value + 1d;
            var b = value - 1d;
            a = Complex.Sqrt(0.5 * a);
            b = Complex.Sqrt(0.5 * b);
            return 2 * Complex.Log(a + b);
        }
        /// <summary>
        /// Returns the angle whose hyperbolic sine is the specified <see cref="Complex"/> number.
        /// </summary>
        /// <param name="value">A <see cref="Complex"/> number representing a hyporbolic sine.</param>
        /// <returns>An angle, mesured in radians, which is the arc hyperbolic sine of <paramref name="value"/></returns>
        public static Complex Asinh(Complex value)
        {
            var d = value.Real - value.Imaginary;
            var re = System.Math.FusedMultiplyAdd(d, d, 1d);
            var im = 2 * value.Real * value.Imaginary;
            Complex t = new(re, im);
            t = Complex.Sqrt(t);
            return Complex.Log(t + value);
        }
        /// <summary>
        /// Returns the angle whose hyperbolic tangens is the specified <see cref="Complex"/> number.
        /// </summary>
        /// <param name="value">A <see cref="Complex"/> number representing a hyporbolic tangens.</param>
        /// <returns>An angle, mesured in radians, which is the arc hyperbolic tangens of <paramref name="value"/></returns>
        public static Complex Atanh(Complex value)
        {
            var i2 = value.Imaginary * value.Imaginary;
            var x = 1 - i2 - value.Real * value.Real;

            var num = 1 + value.Real;
            var den = 1 - value.Real;

            num = System.Math.FusedMultiplyAdd(num, num, i2);
            den = System.Math.FusedMultiplyAdd(den, den, i2);

            var re = 0.25 * (System.Math.Log(num) - System.Math.Log(den));
            var im = 0.5 * System.Math.Atan2(2 * value.Imaginary, x);
            return new(re, im);
        }
        /// <summary>
        /// The norm calculate by this function is also known as <see href="https://en.wikipedia.org/wiki/Field_norm">form norm</see> or <see href="http://mathworld.wolfram.com/AbsoluteSquare.html">absolute square</see>
        /// </summary>
        /// <param name="value">A <see cref="Complex"/> number</param>
        /// <returns>Returns the squared magnitude of the <see cref="Complex"/> number <paramref name="value"/></returns>
        public static double Norm(Complex value)
        {
            var abs = Complex.Abs(value);
            return abs * abs;
        }
        /// <summary>
        /// Returns the specified base logarithm of a specified complex number.
        /// </summary>
        /// <param name="value">A value of logarithm</param>
        /// <param name="base">A base of logarithm</param>
        /// <returns>The <paramref name="base"/> logarithm of <paramref name="value"/> number.</returns>
        public static Complex Log(Complex value, Complex @base) => Complex.Log(value) / Complex.Log(@base);
        /// <summary>
        /// Returns the binary logarithm of a specified complex number.
        /// </summary>
        /// <param name="value">A number</param>
        /// <returns>The binary logarithm of <paramref name="value"/> number</returns>
        public static Complex Log2(Complex value) => Complex.Log(value) / LN2;
        /// <summary>
        /// Returns the projection of the <see cref="Complex"/> number <paramref name="value"/> onto the <see href="https://en.wikipedia.org/wiki/Riemann_sphere">Riemann sphere</see>
        /// </summary>
        /// <param name="value">A <see cref="Complex"/> value</param>
        /// <returns>The projection of <paramref name="value"/> onto the Riemann sphere</returns>
        public static Complex Projection(Complex value)
        {
            var den = Norm(value) + 1;
            var re = 2 * value.Real / den;
            var im = 2 * value.Imaginary / den;
            return new(re, im);
        }
        /// <summary>
        /// Rotate <see cref="Complex"/> number <paramref name="value"/> counter-clockwise by angle <paramref name="radians"/>.
        /// </summary>
        /// <param name="value">A number</param>
        /// <param name="radians">An angle in radians</param>
        /// <returns>Rotated <see cref="Complex"/> number.</returns>
        public static Complex Rotate(Complex value, double radians)
        {
            var r = Complex.Abs(value);
            var arg = value.Phase;
            return Complex.FromPolarCoordinates(r, arg + radians);
        }
        /// <summary>
        /// The signum function generalized to <see cref="Complex"/> numbers.
        /// </summary>
        /// <param name="value">A number</param>
        /// <returns>The signum of <paramref name="value"/></returns>
        public static Complex Sign(Complex value) => (value.Real == 0 && value.Imaginary == 0) ? Complex.Zero : value / value.Magnitude;
        /// <summary>
        /// The signum function generalized to <see cref="Complex"/> numbers.
        /// </summary>
        /// <param name="value">A signed number.</param>
        /// <returns></returns>
        /// <exception cref="ArithmeticException"><inheritdoc cref="System.Math.Sign(double)"/></exception>
        public static int ComplexSign(Complex value)
        {
            var s = System.Math.Sign(value.Real);
            if (s == 0)
            {
                return System.Math.Sign(value.Imaginary);
            }
            else
            {
                return s;
            }
        }
        /// <summary>
        /// <see cref="System.Math.Floor(double)">System.Math.Floor</see> function generalization for <see cref="Complex"/>. Implemented the same way as <see href="https://mathworld.wolfram.com/FloorFunction.html">Floor</see> in the <see href="http://www.wolfram.com/language/">Wolfram Language</see>
        /// </summary>
        /// <param name="value">A <see cref="Complex"/> number</param>
        /// <returns>A new <see cref="Complex"/> object where <br/><see cref="Complex.Real">Real</see> equals <see cref="System.Math.Floor(double)">System.Math.Floor</see>(<paramref name="value"/>.<see cref="Complex.Real">Real</see>) and <br/><see cref="Complex.Imaginary">Imaginary</see> equals <see cref="System.Math.Floor(double)">System.Math.Floor</see>(<paramref name="value"/>.<see cref="Complex.Imaginary">Imaginary</see>)</returns>
        public static Complex Floor(Complex value) => new(System.Math.Floor(value.Real), System.Math.Floor(value.Imaginary));
        /// <summary>
        /// <see cref="System.Math.Ceiling(double)">System.Math.Ceiling</see> function generalization for <see cref="Complex"/>. Implemented the same way as <see href="https://mathworld.wolfram.com/CeilingFunction.html">Ceiling</see> in the <see href="http://www.wolfram.com/language/">Wolfram Language</see>
        /// </summary>
        /// <param name="value">A <see cref="Complex"/> number</param>
        /// <returns>A new <see cref="Complex"/> object where <br/><see cref="Complex.Real">Real</see> equals <see cref="System.Math.Ceiling(double)">System.Math.Ceiling</see>(<paramref name="value"/>.<see cref="Complex.Real">Real</see>) and <br/><see cref="Complex.Imaginary">Imaginary</see> equals <see cref="System.Math.Ceiling(double)">System.Math.Ceiling</see>(<paramref name="value"/>.<see cref="Complex.Imaginary">Imaginary</see>)</returns>
        public static Complex Ceiling(Complex value) => new(System.Math.Ceiling(value.Real), System.Math.Ceiling(value.Imaginary));
        /// <summary>
        /// <see cref="System.Math.Truncate(double)">System.Math.Truncate</see> function generalization for <see cref="Complex"/>. Implemented the same way as <see href="https://mathworld.wolfram.com/IntegerPart.html">IntegerPart</see> in the <see href="http://www.wolfram.com/language/">Wolfram Language</see>
        /// </summary>
        /// <param name="value">A <see cref="Complex"/> number</param>
        /// <returns>A new <see cref="Complex"/> object where <br/><see cref="Complex.Real">Real</see> equals <see cref="System.Math.Truncate(double)">System.Math.Truncate</see>(<paramref name="value"/>.<see cref="Complex.Real">Real</see>) and <br/><see cref="Complex.Imaginary">Imaginary</see> equals <see cref="System.Math.Truncate(double)">System.Math.Truncate</see>(<paramref name="value"/>.<see cref="Complex.Imaginary">Imaginary</see>)</returns>
        public static Complex Truncate(Complex value) => new(System.Math.Truncate(value.Real), System.Math.Truncate(value.Imaginary));
        /// <summary>
        /// <see cref="System.Math.Round(double)">System.Math.Truncate</see> function generalization for <see cref="Complex"/>. Implemented the same way as <see href="https://mathworld.wolfram.com/IntegerPart.html">IntegerPart</see> in the <see href="http://www.wolfram.com/language/">Wolfram Language</see>
        /// </summary>
        /// <param name="value">A <see cref="Complex"/> number</param>
        /// <returns>A new <see cref="Complex"/> object where <br/><see cref="Complex.Real">Real</see> equals <see cref="System.Math.Round(double)">System.Math.Round</see>(<paramref name="value"/>.<see cref="Complex.Real">Real</see>) and <br/><see cref="Complex.Imaginary">Imaginary</see> equals <see cref="System.Math.Round(double)">System.Math.Round</see>(<paramref name="value"/>.<see cref="Complex.Imaginary">Imaginary</see>)</returns>
        public static Complex Round(Complex value) => new(System.Math.Round(value.Real), System.Math.Round(value.Imaginary));
        /// <summary>
        /// Calculate factorial of number <paramref name="n"/>
        /// </summary>
        /// <param name="n">Number of which factorial should be calculated</param>
        /// <returns>Value of <paramref name="n"/>! as <see cref="BigInteger">BigInteger</see></returns>
        /// <exception cref="ArgumentOutOfRangeException">When <paramref name="n"/> is negative</exception>
        public static BigInteger Factorial(int n)
        {
            if (n < 0)
                throw new ArgumentOutOfRangeException(nameof(n), "Cannot be negative");
            else if (n < 1)
                return n;
            else
            {
                BigInteger fact = 1;
                for (int i = 1; i <= n; ++i)
                    fact *= i;
                return fact;
            }
        }
        /// <summary>
        /// Calculate number to roman numerical
        /// </summary>
        /// <param name="number">Number to calculate</param>
        /// <returns>Roman number</returns>
        [CLSCompliant(false)]
        public static string ToRoman(uint number)
        {
            string value = string.Empty;
            for(int i = 0; i < nums.Length && number != 0; ++i)
                while(number >= nums[i])
                {
                    number -= nums[i];
                    value += rum[i];
                }
            return value;
        }
        /// <summary>
        /// Calculate Roman number to normal
        /// </summary>
        /// <param name="roman">Roman number</param>
        /// <returns>uint equivalent</returns>
        [CLSCompliant(false)]
        public static uint ParseRoman(string roman)
        {
            /* Make the input string upper-case,
             * because the dictionary doesn't support lower-case characters. */
            roman = roman.ToUpper();
            /* total = the current total value that will be returned.
             * minus = value to subtract from next numeral. */
            uint total = 0, minus = 0;
            for(int i = 0; i < roman.Length; ++i)
            {
                // Get the value for the current numeral. Takes subtraction into account.
                uint this_numeral = roman_dict[roman[i]] - minus;
                /* Checks if this is the last character in the string, or if the current numeral
                 * is greater than or equal to the next numeral. If so, we will reset our minus
                 * variable and add the current numeral to the total value. Otherwise, we will
                 * subtract the current numeral from the next numeral, and continue. */
                if (i >= roman.Length - 1 || this_numeral + minus >= roman_dict[roman[i + 1]])
                {
                    total += this_numeral;
                    minus = 0;
                }
                else
                {
                    minus = this_numeral;
                }
            }
            return total;
        }
        /// <summary>
        /// Calculate inverse square root of <paramref name="x"/>
        /// </summary>
        /// <param name="x">Non-negative integer</param>
        /// <returns>Integer square root of <paramref name="x"/></returns>
        public static BigInteger InverseSqrt(BigInteger x)
        {
            BigInteger q = 1, r = 0, t;
            while (q <= x)
                q <<= 2;
            while (q > 1)
            {
                q >>= 2;
                t = x - r - q;
                r >>= 1;
                if (t >= 0)
                {
                    x = t;
                    r += q;
                }
            }
            return r;
        }
        /// <summary>
        /// Computes Jacobi symbol of <paramref name="n"/> and <paramref name="k"/>. <paramref name="k"/> should be odd and must be positive
        /// </summary>
        /// <param name="n">Value of n</param>
        /// <param name="k">Value of k. Should be odd and must be positive.</param>
        /// <returns>Jacobi symbol</returns>
        /// <exception cref="ArgumentException"></exception>
        public static BigInteger Jacobi(BigInteger n, BigInteger k)
        {
            if (k < 0 || ((k & 1) == 1)) throw new ArgumentException(null, nameof(k));
            n %= k;
            BigInteger t = 1;
            while(!n.IsZero)
            {
                while((n & 1) == 0)
                {
                    n /= 2;
                    var r = k % 8;
                    if (r == 3 || r == 5)
                        t = -t;
                }
                (k, n) = (n, k);
                if ((n & 4) == 3 && ((k & 3) == 3))
                    t = -t;
                n %= k;
            }
            return k.IsOne ? t : BigInteger.Zero;
        }
        /// <summary>
        /// Computes Jacobi symbol of <paramref name="n"/> and <paramref name="k"/>. <paramref name="k"/> should be odd and must be positive
        /// </summary>
        /// <param name="n">Value of n</param>
        /// <param name="k">Value of k. Should be odd and must be positive.</param>
        /// <returns>Jacobi symbol</returns>
        /// <exception cref="ArgumentException"></exception>
        public static long Jacobi(long n, long k)
        {
            if (k < 0 || ((k & 1) == 1)) throw new ArgumentException(null, nameof(k));
            n %= k;
            long t = 1;
            while (n != 0)
            {
                while ((n & 1) == 0)
                {
                    n /= 2;
                    var r = k % 8;
                    if (r == 3 || r == 5)
                        t = -t;
                }
                (k, n) = (n, k);
                if ((n & 4) == 3 && ((k & 3) == 3))
                    t = -t;
                n %= k;
            }
            return k == 1 ? t : 0;
        }
        /// <summary>
        /// Calculate divisors of given number <paramref name="dividend"/>
        /// </summary>
        /// <param name="dividend"></param>
        /// <returns>Dictionary of number divisors with quantity</returns>
        public static Dictionary<int, int> PrimeDivisors(int dividend)
        {
            Dictionary<int, int> result = new();
            if (dividend < 0) dividend = -dividend;

            if (dividend != 0)
            {
                while(dividend != 1)
                {
                    for(int i = 2; i < dividend; ++i)
                    {
                        if((dividend % i) == 0)
                        {
                            result[i] = result.ContainsKey(i) ? result[i] + 1 : 1;
                            dividend /= i;
                            break;
                        }
                    }
                }
            }
            return result;
        }
        /// <summary>
        /// Calculate divisors of given number <paramref name="dividend"/>
        /// </summary>
        /// <param name="dividend"></param>
        /// <returns>Dictionary of number divisors with quantity</returns>
        public static Dictionary<BigInteger, int> PrimeDivisors(BigInteger dividend)
        {
            Dictionary<BigInteger, int> result = new();
            if (dividend < 0) dividend = -dividend;

            if (dividend != 0)
            {
                while (dividend != 1)
                {
                    for (int i = 2; i < dividend; ++i)
                    {
                        if ((dividend % i) == 0)
                        {
                            result[i] = result.ContainsKey(i) ? result[i] + 1 : 1;
                            dividend /= i;
                            break;
                        }
                    }
                }
            }
            return result;
        }
        /// <summary>
        /// Finds the greatest common divisor of two <see cref="int"/> values.
        /// </summary>
        /// <param name="a">The first value</param>
        /// <param name="b">The second value</param>
        /// <returns>The greatest common division of <paramref name="a"/> and <paramref name="b"/></returns>
        public static int GreatestCommonDivisor(int a, int b)
        {
            while (b != 0) b = a % (a = b);
            return a;
        }
        /// <summary>
        /// Finds the greatest common divisor of two <see cref="long"/> values.
        /// </summary>
        /// <param name="a">The first value</param>
        /// <param name="b">The second value</param>
        /// <returns>The greatest common division of <paramref name="a"/> and <paramref name="b"/></returns>
        public static long GreatestCommonDivisor(long a, long b)
        {
            while (b != 0) b = a % (a = b);
            return a;
        }
        /// <summary>
        /// The Minkowski question-mark function converts the continued fraction representation [a0; a1, a2, a3, ...] of a number into a binary decimal representation in which the integer part a0 is unchanged and the a1, a2, ... become alternating runs of binary zeroes and ones of those lengths. The decimal point takes the place of the first zero.
        /// </summary>
        /// <param name="x">Value to converts</param>
        /// <returns>Converted value</returns>
        public static double Minkowski(double x)
        {
            long p = (long)x;
            if (p > x) --p;
            long q = 1, r = p + 1, s = 1, m, n;
            double d = 1, y = p;
            if (x < p || (p < 0) ^ (r <= 0))
                return x;
            for(; ; )/* invariants: q * r - p * s == 1 && (double)p / q <= x && x < (double)r / s */
            {
                d /= 2;
                if (y + d == y)
                    break; /* reached max possible precision */
                m = p + r;
                if ((m < 0) ^ (p < 0))
                    break; /* sum overflowed */
                n = q + s;
                if (n < 0)
                    break; /* sum overflowed */

                if (x < (double)m / n)
                {
                    r = m;
                    s = n;
                }
                else
                {
                    y += d;
                    p = m;
                    q = n;
                }
            }
            return y + d; /* final round-off */
        }
    }
}
