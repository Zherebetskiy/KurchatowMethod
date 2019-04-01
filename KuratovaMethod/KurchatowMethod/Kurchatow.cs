using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using System;

namespace KurchatowMethod
{
    public class Kurchatow
    {
        private readonly bool useMatrix = false;
        private double[,] A;
        private double[] B;
        private readonly int n;
        private readonly double Eps;

        public int IterationCount { get; set; }

        public double[] X0 { get; set; }
        public double[] X1 { get; set; }

        public Kurchatow(int n, double Eps)
        {
            this.Eps = Eps;
            this.n = n;
            A = new double[n, n];
            B = new double[n];
            X0 = new double[n];
            X1 = new double[n];
            useMatrix = n > 1;
        }

        public double Func(double x)
        {
            return Math.Pow(x, 3) - 15 * x + 10;
            //return Math.Pow(x, 3) - 2 * Math.Sin(x) - 1;
        }

        public double FuncPoh(double x)
        {
            return 3 * Math.Pow(x, 2) - 15;
            //return Math.Pow(x, 2) + 2 * Math.Cos(x);
        }

        public double Sys1(double x, double y)
        {
            return Math.Cos(x - 1) + y - 1;
            //return 2 * x - Math.Sin(0.5 * (x - y));
        }

        public double Sys2(double x, double y)
        {
            return Math.Sin(x - 1) + y - 1;
            //return 2 * y - Math.Cos(0.5 * (x + y));
        }

        private void MatrixBuilder(double[] x_prev, double[] x_next)
        {
            double[] x = new double[n];
            for (int i = 0; i < n; i++)
            {
                x[i] = 2 * x_next[i] - x_prev[i];
            }

            if (useMatrix)
            {
                A[0, 0] = 1.0 / (x[0] - x_prev[0]) * (Sys1(x[0], x_prev[1]) - Sys1(x_prev[0], x_prev[1]));
                A[0, 1] = 1.0 / (x[1] - x_prev[1]) * (Sys1(x[0], x[1]) - Sys1(x[0], x_prev[1]));
                A[1, 0] = 1.0 / (x[0] - x_prev[0]) * (Sys2(x[0], x_prev[1]) - Sys2(x_prev[0], x_prev[1]));
                A[1, 1] = 1.0 / (x[1] - x_prev[1]) * (Sys2(x[0], x[1]) - Sys2(x[0], x_prev[1]));
            }
            else
            {
                A[0, 0] = 1.0 / (x[0] - x_prev[0]) * (Func(x[0]) - Func(x_prev[0]));
            }
        }

        private void VectorBuilder(double[] x)
        {
            if (useMatrix)
            {
                B[0] = -Sys1(x[0], x[1]);
                B[1] = -Sys2(x[0], x[1]);
            }
            else
            {
                B[0] = -Func(x[0]);
            }
        }

        public double[] Iteration(double[] x, double[] x_next)
        {
            IterationCount = 0;
            double[] x_prev = new double[n];
            double[] sigma = new double[n];

            do
            {
                x_prev = ChangeVector(x);
                x = ChangeVector(x_next);
                sigma = GetSigma(x_prev, x);

                for (int i = 0; i < n; i++)
                {
                    x_next[i] = x[i] + sigma[i];
                }
                IterationCount++;
            }
            while (Norma(x, x_next) > Eps);

            return x_next;
        }

        public void InitialPraram(double[] x0, double jot)
        {
            if (useMatrix)
            {
                for (int i = 0; i < n; i++)
                {
                    X0[i] = x0[i];
                    X1[i] = X0[i] - 0.0000001;
                }
            }
            else
            {
                var f = Func(x0[0]);
                var fp = FuncPoh(x0[0]);
                var beta = 1.0 / Math.Abs(FuncPoh(x0[0]));
                var nue = Math.Abs(Func(x0[0]) / FuncPoh(x0[0]));
                var alpha = beta * nue * jot;
                double r = 0.0;

                if (alpha < 0.5)
                {
                    r = (1 + Math.Sqrt(1 - 2 * alpha)) / (beta * jot);
                }
                else
                {
                    r = (1 - Math.Sqrt(1 - 2 * alpha)) / (beta * jot);
                }

                var a = x0[0] - r;
                var b = x0[0] + r;


                for (int i = 0; i < n; i++)
                {
                    X0[i] = x0[i];
                    X1[i] = X0[i] - 0.0000001;
                }

                Console.WriteLine($"a = {a} b = {b}");

                for (int i = 0; i < n; i++)
                {
                    Console.WriteLine($"X0[{i}] = {X0[i]} X1[{i}]  = {X1[i]}");
                }

                Console.WriteLine($"Beta = {beta} Nue = {nue} Jot = {jot}");
                Console.WriteLine($"Alpha = {alpha}");
                Console.WriteLine($"R = {r}");
            }
        }

        private double[] ChangeVector(double[] x)
        {
            double[] y = new double[n];

            for (int i = 0; i < n; i++)
            {
                y[i] = x[i];
            }

            return y;
        }

        private double[] GetSigma(double[] x_prev, double[] x_next)
        {
            MatrixBuilder(x_prev, x_next);
            VectorBuilder(x_next);

            Matrix<double> matrix = DenseMatrix.OfArray(A);
            Vector<double> vector = DenseVector.OfArray(B);

            return matrix.Solve(vector).AsArray();
        }

        private double Norma(double[] x1, double[] x2)
        {
            if (useMatrix)
            {
                return Math.Max(Math.Abs(x1[0] - x2[0]), Math.Abs(x1[1] - x2[1]));
            }

            return Math.Abs(x1[0] - x2[0]);
        }
    }
}
