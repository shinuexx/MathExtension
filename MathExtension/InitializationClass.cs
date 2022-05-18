using System.Reflection;
using System.Runtime.InteropServices;

namespace ShInUeXx.Numerics
{
    internal static class InitializationClass
    {
        internal static bool IsInitialized { get; private set; } = false;
        internal static T? Initialize<T>() where T : struct
        {
            if (!IsInitialized)
            {
                IsInitialized = true;
                NativeLibrary.SetDllImportResolver(Assembly.GetExecutingAssembly(), DllImportResolver);
            }
            return null;
        }
        internal static IntPtr DllImportResolver(string library_name, Assembly assembly, DllImportSearchPath? search_path)
        {
            if (library_name == "cpp-math")
            {
                if (RuntimeInformation.IsOSPlatform(OSPlatform.Windows))
                {
                    return NativeLibrary.Load("ucrtbase.dll", assembly, search_path);
                }
                else if (RuntimeInformation.IsOSPlatform(OSPlatform.Linux) || RuntimeInformation.IsOSPlatform(OSPlatform.FreeBSD))
                {
                    return NativeLibrary.Load("libm.so.6", assembly, search_path);
                }
                else if (RuntimeInformation.IsOSPlatform(OSPlatform.OSX))
                {
                    return NativeLibrary.Load("libSystem.B.dylib", assembly, search_path);
                }
            }
            return IntPtr.Zero;
        }
    }
}