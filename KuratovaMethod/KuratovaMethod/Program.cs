using System;
using KurchatowMethod;

namespace KuratovaMethod
{
    class Program
    {
        static void Main(string[] args)
        {
            var n = 1;
            var x0 = new double[] { 4,2 };
            var jot = 27;
            var eps = 0.000001;
            Kurchatow kurchatow = new Kurchatow(n, eps);
            kurchatow.InitialPraram(x0, jot);
            //var res = kurchatow.Iteration(new double[] {0, 0.5 }, new double[] { -0.16, 0.49});
            //var res = kurchatow.Iteration(new double[] { 2 }, new double[] { 1.59 });
            var res = kurchatow.Iteration(new double[] { kurchatow.X0[0] }, new double[] { kurchatow.X1[0]});

            Console.WriteLine($"Iteration = {kurchatow.IterationCount}");
            for (int i = 0; i < n; i++)
            {
                Console.WriteLine($"res[{i}] = {res[i]}");
            }

            if (n > 1)
            {
                Console.WriteLine($" Func1 = {kurchatow.Sys1(res[0], res[1])}  Func2 = {kurchatow.Sys2(res[0], res[1])}");
            }
            else
            {
                Console.WriteLine($" Func(res) = {kurchatow.Func(res[0])}");
            }


            Console.Read();
        }
    }
}
