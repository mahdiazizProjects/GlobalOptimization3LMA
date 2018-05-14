using System;
using System.Collections;
public class function_evaluation2013
{
    private double ub;
    private double lb;
    public static double[,] Eye(int N)
    {
        double[,] results = new double[N, N];
        for(int ii=0;ii<N;ii++)
            for(int jj=0;jj<N;jj++)
            {
                if (ii == jj)
                    results[ii, jj] = 1;
            }
        return results;
    }
    public double[,] MultiplyMatrix(double[,] a, double[,] b)
    {
        double[,] c = new double[a.GetLength(0), b.GetLength(1)];
        if (a.GetLength(1) == b.GetLength(0))
        {
            for (int i = 0; i < c.GetLength(0); i++)
            {
                for (int j = 0; j < c.GetLength(1); j++)
                {
                    c[i, j] = 0;
                    for (int k = 0; k < a.GetLength(1); k++) // OR k<b.GetLength(0)
                        c[i, j] = c[i, j] + a[i, k] * b[k, j];
                }
            }
        }
        return c;
    }
    public double replace(double[]m)
    {
        double sum=0;
        for (int ii = 0; ii < m.Length; ii++)
            sum+=m[ii] * m[ii];
        return sum;
    }
    public double[,]replace2d(double[,]m)
    {
        double[,]result=new double[m.GetLength(1),m.GetLength(0)];
        for (int ii = 0; ii < m.GetLength(0); ii++)
            for (int jj = 0; jj < m.GetLength(1); jj++)
                result[jj, ii] = m[ii, jj];
        return result;
    }
    public double Update_evolution_paths(double cs,ref double[] ps, double mueff, double[,] invsqrtC, double[,] xmean, double[,] xold, double sigma, int counteval, double lambda, double chiN,int N,double cc,ref double[]pc)
    {
        double[,] diff = new double[xmean.GetLength(0), 1];
        double coeff = Math.Sqrt(cs * (2 - cs) * mueff);
        for (int ii = 0; ii < diff.GetLength(0); ii++)
        {
            diff[ii, 0] = xmean[ii, 0] - xold[ii, 0];
        }
        double[,] invmuldiff = MultiplyMatrix(invsqrtC, diff);
        double[,] invmuldiffcoef = new double[invmuldiff.GetLength(0), invmuldiff.GetLength(1)];
        for (int ii = 0; ii < invmuldiff.GetLength(0); ii++)
            for (int jj = 0; jj < invmuldiff.GetLength(1); jj++)
                invmuldiffcoef[ii, jj] = coeff * invmuldiff[ii, jj] / sigma;
        for (int ii = 0; ii < ps.Length; ii++)
        {
            ps[ii] = (1 - cs) * ps[ii] + invmuldiffcoef[ii, 0];
        }
        ////////////////norm ps
        double norm=0;
        for(int ii=0;ii<ps.Length;ii++)
        {
            norm+=Math.Pow(ps[ii] , 2);
        }
        norm = Math.Sqrt(norm);
        double hsig=-1;
        if((norm/Math.Sqrt(Math.Pow(1-(1-cs),(2*counteval/lambda)))/chiN) < (1.4 + 2/(N+1)))
            hsig =1;
        else
            hsig=0;
        for (int ii = 0; ii < pc.Length; ii++)
            pc[ii] = (1 - cc) * pc[ii] + hsig * Math.Sqrt(cc * (2 - cc) * mueff) * diff[ii, 0] / sigma;
        return hsig;
    }
    public double[,]diags(double[,] weights)
    {
        double[,] result = new double[weights.GetLength(0), weights.GetLength(0)];
        for (int ii = 0; ii < weights.Length; ii++)
            for (int jj = 0; jj < weights.Length; jj++)
                if (ii == jj)
                    result[ii, jj] = weights[ii,0];
            
        return result;
    }
    public double[] create_new_gen_cma_es(int N, double[,] xmean,double sigma,double[,]B,double[,]D)
    {
        // arx(:,k) = xmean + sigma * B * (D .* randn(N,1)); % m + sig * Normal(0,C) 
        double[] result = new double[N];
        Random rr=new Random();
        double[,]Dr=new double[N,1];
        for (int ii = 0; ii < D.Length;ii++ )
        {
            Dr[ii,0] = D[ii,0] * Normal_Distribution(rr);
        }
        double[,]mulbdr=MultiplyMatrix(B, Dr);
        for (int ii = 0; ii < xmean.GetLength(0);ii++ )
        {
            result[ii] = xmean[ii,0] * sigma * mulbdr[ii, 0];
        }
        return result;
    }
    public function_evaluation2013(double ub, double lb, int D)
    {
        this.ub = ub;
        this.lb = lb;
    }

    public function_evaluation2013()
    {
        // TODO: Complete member initialization
    }
    private double f1(double[] x, int D, double[] xopt)
    {
        double[] y = new double[x.Length];
        for (int ii = 0; ii < x.Length; ii++)
        {
            y[ii] = x[ii] - xopt[ii];
        }
        return elliptic(y, D);
    }
    public double Normal_Distribution(Random rand)
    {
        double u1 = rand.NextDouble(); //these are uniform(0,1) random doubles
        double u2 = rand.NextDouble();
        double randStdNormal = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Sin(2.0 * Math.PI * u2); //random normal(0,1)
        return randStdNormal;
    }
    //private double elliptic(double[]x,int D,double condition)
    //{
    //    double[]lin=Generate.LinearSpaced(D,0, 1);
    //    double[]coefficients=new double[lin.Length];
    //    for(int ii=0;ii<lin.Length;ii++)
    //    {
    //       coefficients[ii] =  Math.Pow(condition, lin[ii]);
    //    }
    //    double[] xtemp = T_irreg(x);
    //    double fit = 0;
    //    for(int ii=0;ii<x.Length;ii++)
    //    {
    //      fit+= Math.Pow(coefficients[ii] * xtemp[ii],2);
    //    }
    //    return fit;
    //}
    //private double[] T_irreg(double []f)
    //{
    //    double a = 0.1;
    //    double[]g =(double[])f.Clone(); 
    //    ArrayList idx=new ArrayList();
    //    for(int uu=0;uu<f.Length;uu++)
    //    {
    //        if(f[uu]>0)
    //            idx.Add(uu);
    //    }
    //    for(int ii=0;ii<idx.Count;ii++)
    //    {
    //        g[(int)idx[ii]] = Math.Log(f[(int)idx[ii]]) / a;
    //        g[(int)idx[ii]] = Math.Pow(Math.Exp(g[(int)idx[ii]] + 0.49 * (Math.Sin(g[(int)idx[ii]]) + Math.Sin(0.79 * g[(int)idx[ii]]))), a);
    //    }
    //    idx.Clear();
    //    for(int uu=0;uu<f.Length;uu++)
    //    {
    //        if(f[uu]<0)
    //            idx.Add(uu);
    //    }
    //    for(int ii=0;ii<idx.Count;ii++)
    //    {
    //        g[(int)idx[ii]] = Math.Log(-f[(int)idx[ii]]) / a;
    //        g[(int)idx[ii]] = -Math.Pow(Math.Exp(g[(int)idx[ii]] + 0.49 * (Math.Sin(g[(int)idx[ii]]) + Math.Sin(0.79 * g[(int)idx[ii]]))), a);
    //    }
    //    return g;
    //}
    private double elliptic(double[] x, int D)
    {
        double result = 0.0;
        int i;

        transform_osz(ref x);

        // for(i = dim - 1; i >= 0; i--) {
        for (i = 0; i < D; i++)
        {
            // printf("%f\n", pow(1.0e6,  i/((double)(dim - 1)) ));
            result += Math.Pow(1.0e6, i / ((double)(D - 1))) * x[i] * x[i];
        }

        return (result);
    }
    private double elliptic(double[] x, double dim)
    {
        double result = 0.0;
        int i;

        transform_osz(ref x);

        // for(i = dim - 1; i >= 0; i--) {
        for (i = 0; i < dim; i++)
        {
            // printf("%f\n", pow(1.0e6,  i/((double)(dim - 1)) ));
            result += Math.Pow(1.0e6, i / ((double)(dim - 1))) * x[i] * x[i];
        }

        return (result);
    }
    private void transform_osz(ref double[] f)
    {
        // apply osz transformation to z
        for (int i = 0; i < f.Length; ++i)
        {
            f[i] = Math.Sign(f[i]) * Math.Exp(hat(f[i]) + 0.049 * (Math.Sin(c1(f[i]) * hat(f[i])) + Math.Sin(c2(f[i]) * hat(f[i]))));
        }
    }
    private double f2(double[] x, int D, double[] xopt)
    {

        double result;
        double[] y = new double[x.Length];

        for (int i = 0; i < D; i++)
        {
            y[i] = x[i] - xopt[i];
        }

        // // T_{osz}
        // transform_osz(anotherz, dimension);

        // // T_{asy}^{0.2}
        // transform_asy(anotherz, 0.2);

        // // lambda
        // Lambda(anotherz, 10);

        result = rastrigin(y, D);
        return (result);
    }

    private double rastrigin(double[] x, int D)
    {
        double sum = 0;
        // T_{osz}
        transform_osz(ref x);

        // T_{asy}^{0.2}
        transform_asy(ref x, 0.2, D);

        // lambda
        Lambda(ref x, 10, D);

        for (int i = D - 1; i >= 0; i--)
        {
            sum += x[i] * x[i] - 10.0 * Math.Cos(2 * Math.PI * x[i]) + 10.0;
        }

        return (sum);
    }
    /*------------------------------------------------------------------------------
% This transformation function is used to break the symmetry of symmetric 
% functions.
%------------------------------------------------------------------------------*/
    private double[] T_asy(double[] f, double beta, int D)
    {
        double[] g = (double[])f.Clone();
        for (int i = 0; i < D; ++i)
        {
            if (f[i] > 0)
            {
                g[i] = Math.Pow(f[i], 1 + beta * i / ((double)(D - 1)) * Math.Sqrt(f[i]));
            }
        }
        return g;
    }
    private double hat(double x)
    {
        if (x == 0)
        {
            return 0;
        }
        else
        {
            return Math.Log(Math.Abs(x));
        }
    }
    private double c1(double x)
    {
        if (x > 0)
        {
            return 10;
        }
        else
        {
            return 5.5;
        }
    }

    private double c2(double x)
    {
        if (x > 0)
        {
            return 7.9;
        }
        else
        {
            return 3.1;
        }
    }
    private void transform_asy(ref double[] f, double beta, int D)
    {
        for (int i = 0; i < D; ++i)
        {
            if (f[i] > 0)
            {
                f[i] = Math.Pow(f[i], 1 + beta * i / ((double)(D - 1)) * Math.Sqrt(f[i]));
            }
        }
    }
    private void Lambda(ref double[] z, double alpha, int D)
    {
        for (int i = 0; i < D; ++i)
        {
            z[i] = z[i] * Math.Pow(alpha, 0.5 * i / ((double)(D - 1)));
        }
    }
    private double F3(double[] x, double[] xopt, int D)
    {
        int i;
        double result;
        double[] y = new double[x.Length];
        for (i = D - 1; i >= 0; i--)
        {
            y[i] = x[i] - xopt[i];
        }

        result = ackley(y, D);

        return (result);
    }
    private double ackley(double[] x, int D)
    {
        double sum1 = 0.0;
        double sum2 = 0.0;
        double sum;
        int i;

        // T_{osz}
        transform_osz(ref x);

        // T_{asy}^{0.2}
        transform_asy(ref x, 0.2, D);

        // lambda
        Lambda(ref x, 10, D);

        for (i = D - 1; i >= 0; i--)
        {
            sum1 += (x[i] * x[i]);
            sum2 += Math.Cos(2.0 * Math.PI * x[i]);
        }

        sum = -20.0 * Math.Exp(-0.2 * Math.Sqrt(sum1 / D)) - Math.Exp(sum2 / D) + 20.0 + Math.E;
        return (sum);
    }
    private double F4(double[] x, double[] xopt, int[] perm, double[] r25, double[] r50, double[] r100, int[] s, double[] w, int D)
    {
        int i;
        double result = 0.0;
        double[] y = new double[x.Length];
        double[] y2 = new double[x.Length];

        for (i = 0; i < D; i++)
        {
            y[i] = x[i] - xopt[i];
        }

        // for (int i = 0; i < dimension; ++i)
        //   {
        //     cout<<anotherz[i]<<endl;
        //   }
        // cout<<endl;

        // // T_{osz}
        // transform_osz(anotherz);

        // s_size non-separable part with rotation
        int c = 0;
        for (i = 0; i < s.Length; i++)
        {
            // cout<<"c="<<c<<", i="<<i<<endl;
            y2 = rotateVector(s, i,ref c, perm, y, y2, r25, r50, r100);
            // cout<<"done rot"<<endl;
            result += w[i] * elliptic(y2, s[i]);
            // cout<<result<<endl;
        }

        // one separable part without rotation
        double[] z = new double[D - c];
        for (i = c; i < D; i++)
        {
            // cout<<i-c<<" "<<Pvector[i]<<" "<<anotherz[Pvector[i]]<<endl;
            z[i - c] = y2[perm[i]-1];
        }

        // cout<<"sep\n"<<elliptic(z, dimension-c)<<endl;

        result += elliptic(z, D - c);

        return (result);
    }
    private double[] rotateVector(int[] s, int i,ref int c, int[] perm, double[] y, double[] y2, double[] r25, double[] r50, double[] r100)
    {
        // cout<<"s "<<s[i]<<endl;

        // copy values into the new vector
        for (int j = c; j < c + s[i]; ++j)
        {
            // cout<<"j-c "<<j-c<<" p "<<Pvector[j]<<endl;
            y2[j - c] = y[perm[j]-1];
        }
        // cout<<"copy done"<<endl;

        if (s[i] == 25)
        {
            y2 = multiply(y2, r25, s[i]);
        }
        else if (s[i] == 50)
        {
            y2 = multiply(y2, r50, s[i]);
        }
        else if (s[i] == 100)
        {
            y2 = multiply(y2, r100, s[i]);
        }
        else
        {

        }
        c = c + s[i];
        return y2;
    }
    private double[] multiply(double[] vector, double[] matrix, int dim)
    {
        int i, j;
        //double*result = (double*)malloc(sizeof(double) * dim);
        double[] result = new double[dim];

        for (i = dim - 1; i >= 0; i--)
        {
            result[i] = 0;

            for (j = dim - 1; j >= 0; j--)
            {
                result[i] += vector[j] * matrix[dim * j + i];
            }
        }

        return (result);
    }
    private double F5(double[] x, double[] xopt, int[] perm, double[] r25, double[] r50, double[] r100, int[] s, double[] w, int D)
    {
        int i;
        double result = 0.0;
        double[] y = new double[x.Length];
        double[] y2 = new double[x.Length];
        for (i = 0; i < D; i++)
        {
            y[i] = x[i] - xopt[i];
        }

        // put them inside rastrigin function
        // // T_{osz}
        // transform_osz(anotherz, dimension);

        // // T_{asy}^{0.2}
        // transform_asy(anotherz, 0.2);

        // // lambda
        // Lambda(anotherz, 10);

        // s_size non-separable part with rotation
        int c = 0;
        for (i = 0; i < s.Length; i++)
        {
            // cout<<"c="<<c<<", i="<<i<<endl;
            y2 = rotateVector(s, i,ref c, perm, y, y2, r25, r50, r100);
            // cout<<"done rot"<<endl;
            result += w[i] * rastrigin(y2, s[i]);
            // cout<<result<<endl;
        }

        // one separable part without rotation
        double[] z = new double[D - c];
        for (i = c; i < D; i++)
        {
            // cout<<i-c<<" "<<Pvector[i]<<" "<<anotherz[Pvector[i]]<<endl;
            z[i - c] = y[perm[i]-1];
        }

        // cout<<"sep"<<endl;
        // cout<<rastrigin(z, dimension-c)<<endl;

        result += rastrigin(z, D - c);
        // free(z);
        return (result);
    }
    private double F6(double[] x, double[] xopt, int[] perm, double[] r25, double[] r50, double[] r100, int[] s, double[] w, int D)
    {
        int i;
        double result = 0.0;
        double[] y = new double[x.Length];
        double[] y2 = new double[x.Length];
        for (i = 0; i < D; i++)
        {
            y[i] = x[i] - xopt[i];
        }

        // cout<<"non"<<endl;

        // s_size non-separable part with rotation
        int c = 0;
        for (i = 0; i < s.Length; i++)
        {
            // cout<<"c="<<c<<", i="<<i<<endl;
            y2 = rotateVector(s, i,ref c, perm, y, y2, r25, r50, r100);
            // cout<<"done rot"<<endl;
            result += w[i] * ackley(y2, s[i]);
            // cout<<result<<endl;
        }

        // one separable part without rotation
        double[] z = new double[D - c];
        for (i = c; i < D; i++)
        {
            // cout<<i-c<<" "<<Pvector[i]<<" "<<anotherz[Pvector[i]]<<endl;
            z[i - c] = y[perm[i]-1];
        }
        result += ackley(z, D - c);
        return (result);
    }
    private double F7(double[] x, double[] xopt, int[] perm, double[] r25, double[] r50, double[] r100, int[] s, double[] w, int D)
    {
        int i;
        double result = 0.0;
        double[] y = new double[x.Length];
        double[] y2 = new double[x.Length];
        for (i = 0; i < D; i++)
        {
            y[i] = x[i] - xopt[i];
        }

        // s_size non-separable part with rotation
        int c = 0;
        for (i = 0; i < s.Length; i++)
        {
            // cout<<"c="<<c<<", i="<<i<<endl;
            y2 = rotateVector(s, i,ref c, perm, y, y2, r25, r50, r100);
            // cout<<"done rot"<<endl;
            result += w[i] * schwefel(y2, s[i]);
            // cout<<result<<endl;
        }

        // one separable part without rotation
        double[] z = new double[D - c];
        for (i = c; i < D; i++)
        {
            // cout<<i-c<<" "<<Pvector[i]<<" "<<anotherz[Pvector[i]]<<endl;
            z[i - c] = y[perm[i]-1];
        }

        result += sphere(z, D - c);
        return (result);
    }
    private double schwefel(double[] x, int dim)
    {
        int j;
        double s1 = 0;
        double s2 = 0;

        // T_{osz}
        transform_osz(ref x);

        // T_{asy}^{0.2}
        transform_asy(ref x, 0.2, dim);

        for (j = 0; j < dim; j++)
        {
            s1 += x[j];
            s2 += (s1 * s1);
        }

        return (s2);
    }
    private double sphere(double[] x, int dim)
    {
        double sum = 0;
        int i;

        for (i = dim - 1; i >= 0; i--)
        {
            sum += Math.Pow(x[i], 2);
        }

        return (sum);
    }
    private double F8(double[] x, double[] xopt, int[] perm, double[] r25, double[] r50, double[] r100, int[] s, double[] w, int D)
    {
        int i;
        double result = 0.0;
        double[] y = new double[x.Length];
        double[] y2 = new double[x.Length];
        for (i = 0; i < D; i++)
        {
            y[i] = x[i] - xopt[i];
        }

        // s_size non-separable part with rotation
        int c = 0;
        for (i = 0; i < s.Length; i++)
        {
            // cout<<"c="<<c<<", i="<<i<<endl;
            y2 = rotateVector(s, i,ref c, perm, y, y2, r25, r50, r100);
            // cout<<"done rot"<<endl;
            result += w[i] * elliptic(y2, s[i]);
            // cout<<result<<endl;
        }

        return (result);
    }
    private double F9(double[] x, double[] xopt, int[] perm, double[] r25, double[] r50, double[] r100, int[] s, double[] w, int D)
    {
        int i;
        double result = 0.0;
        double[] y = new double[x.Length];
        double[] y2 = new double[x.Length];
        for (i = 0; i < D; i++)
        {
            y[i] = x[i] - xopt[i];
        }

        // s_size non-separable part with rotation
        int c = 0;
        for (i = 0; i < s.Length; i++)
        {
            // cout<<"c="<<c<<", i="<<i<<endl;
            y2 = rotateVector(s, i,ref c, perm, y, y2, r25, r50, r100);
            // cout<<"done rot"<<endl;
            result += w[i] * rastrigin(y2, s[i]);
            // cout<<result<<endl;
        }


        return (result);
    }
    private double F10(double[] x, double[] xopt, int[] perm, double[] r25, double[] r50, double[] r100, int[] s, double[] w, int D)
    {
        int i;
        double result = 0.0;
        double[] y = new double[x.Length];
        for (i = 0; i < D; i++)
        {
            y[i] = x[i] - xopt[i];
        }

        // s_size non-separable part with rotation
        int c = 0;
        for (i = 0; i < s.Length; i++)
        {
            // cout<<"c="<<c<<", i="<<i<<endl;
            double[] y2 = new double[s[i]];
            y2 = rotateVector(s, i,ref c, perm, y,y2 ,r25, r50, r100);
            // cout<<"done rot"<<endl;
            result += w[i] * ackley(y2, s[i]);
            // cout<<result<<endl;

        }

        return (result);
    }
    private double F11(double[] x, double[] xopt, int[] perm, double[] r25, double[] r50, double[] r100, int[] s, double[] w, int D)
    {
        int i;
        double result = 0.0;
        double[] y = new double[x.Length];
        double[] y2 = new double[x.Length];
        for (i = 0; i < D; i++)
        {
            y[i] = x[i] - xopt[i];
        }

        // s_size non-separable part with rotation
        int c = 0;
        for (i = 0; i < s.Length; i++)
        {
            // cout<<"c="<<c<<", i="<<i<<endl;
            y2 = rotateVector(s, i,ref c, perm, y, y2, r25, r50, r100);
            // cout<<"done rot"<<endl;
            result += w[i] * schwefel(y2, s[i]);
            // cout<<result<<endl;
        }
        return (result);
    }
    private double F12(double[] x, double[] xopt, int D)
    {
        int i;
        double result = 0.0;

        double[] y = new double[x.Length];
        for (i = 0; i < D; i++)
        {
            y[i] = x[i] - xopt[i];
        }

        result = rosenbrock(y, D);
        return (result);
    }
    private double rosenbrock(double[] x, int dim)
    {
        int j;
        double oz, t;
        double s = 0.0;
        j = dim - 1;

        for (--j; j >= 0; j--)
        {
            oz = x[j + 1];
            t = ((x[j] * x[j]) - oz);
            s += (100.0 * t * t);
            t = (x[j] - 1.0);
            s += (t * t);
        }
        return (s);
    }
    private double F13(double[] x, double[] xopt, int[] perm, double[] r25, double[] r50, double[] r100, int[] s, double[] w, int D, int overlap)
    {
        int i;
        double result = 0.0;
        double[] y = new double[x.Length];
        double[] y2 = new double[x.Length];
        for (i = 0; i < D; i++)
        {
            y[i] = x[i] - xopt[i];
        }

        // s_size non-separable part with rotation
        int c = 0;
        for (i = 0; i < s.Length; i++)
        {
            // cout<<"c="<<c<<", i="<<i<<endl;
            y2 = rotateVectorConform(s, i, c, perm, y, y2, r25, r50, r100, overlap);
            // cout<<"done rot"<<endl;
            result += w[i] * schwefel(y2, s[i]);
            // cout<<result<<endl;
        }

        return (result);
    }
    private double[] rotateVectorConform(int[] s, int i, int c, int[] perm, double[] y, double[] y2, double[] r25, double[] r50, double[] r100, int overlap)
    {
        double[] z = new double[s[i]];
        // printf("i=%d\tc=%d\tl=%d\tu=%d\n", i, c, c - (i)*overlap, c +s[i] - (i)*overlap);

        // copy values into the new vector
        for (int j = c - i * overlap; j < c + s[i] - i * overlap; ++j)
        {
            // cout<<"j-c "<<j-c<<" p "<<Pvector[j]<<endl;
            z[j - (c - i * overlap)] = y[perm[j]-1];
        }
        // cout<<"copy done"<<endl;
        if (s[i] == 25)
        {
            y2 = multiply(z, r25, s[i]);
        }
        else if (s[i] == 50)
        {
            y2 = multiply(z, r50, s[i]);
        }
        else if (s[i] == 100)
        {
            y2 = multiply(z, r100, s[i]);
        }
        else
        {

        }
        c = c + s[i];
        return y2;
    }
    private double F14(double[] x, double[,] xopt, int[] perm, double[] r25, double[] r50, double[] r100, int[] s, double[] w, int D, int overlap)
    {
        int i;
        double result = 0.0;
        double[] y2 = new double[x.Length];

        // s_size non-separable part with rotation
        int c = 0;
        for (i = 0; i < s.Length; i++)
        {

            // cout<<"c="<<c<<", i="<<i<<endl;
            y2 = rotateVectorConflict(s, i, c, perm, x, y2, r25, r50, r100, overlap, xopt);
            // cout<<"done rot"<<endl;
            result += w[i] * schwefel(y2, s[i]);

        }

        return (result);
    }
    private double[] rotateVectorConflict(int[] s, int i, int c, int[] perm, double[] x, double[] y2, double[] r25, double[] r50, double[] r100, int overlap, double[,] xop)
    {
        double[] z = new double[s[i]];
        // printf("i=%d\tc=%d\tl=%d\tu=%d\n", i, c, c - (i)*overlap, c +s[i] - (i)*overlap);

        // copy values into the new vector
        for (int j = c - i * overlap; j < c + s[i] - i * overlap; ++j)
        {
            // cout<<"j-c "<<j-c<<" p "<<Pvector[j]<<endl;
            z[j - (c - i * overlap)] = x[perm[j]] - xop[i, j - (c - i * overlap)];
        }

        // cout<<"copy done"<<endl;
        if (s[i] == 25)
        {
            y2 = multiply(z, r25, s[i]);
        }
        else if (s[i] == 50)
        {
            y2 = multiply(z, r50, s[i]);
        }
        else if (s[i] == 100)
        {
            y2 = multiply(z, r100, s[i]);
        }
        else
        {
        }
        c = c + s[i];
        return y2;
    }
    private double F15(double[] x, double[] xopt, int D)
    {
        int i;
        double result = 0.0;
        double[] y = new double[x.Length];
        for (i = 0; i < D; i++)
        {
            y[i] = x[i] - xopt[i];
        }


        result = schwefel(y, D);

        return (result);
    }
    public double Choosing_Function(int function,double[]x,int D,double[]xopt,int[]perm,double[]r25,double[]r50,double[]r100,int[]s,double[]w,int overlap)
    {
        double result = 0;
        if(function==1)
        {
            result = f1(x, D, xopt);
        }
        else if(function==2)
        {
            result = f2(x, D, xopt);
        }
        else if (function == 3)
        {
            result = F3(x, xopt, D);
        }
        else if (function == 4)
        {
            result = F4(x, xopt, perm, r25, r50, r100, s, w, D);
        }
        else if (function == 5)
        {
            result = F5(x, xopt, perm, r25, r50, r100, s, w, D);
        }
        else if (function == 6)
        {
            result = F6(x, xopt, perm, r25, r50, r100, s, w, D);
        }
        else if (function == 7)
        {
            result = F7(x, xopt, perm, r25, r50, r100, s, w, D);
        }
        else if (function == 8)
        {
            result = F8(x, xopt, perm, r25, r50, r100, s, w, D);
        }
        else if (function == 9)
        {
            result = F9(x, xopt, perm, r25, r50, r100, s, w, D);
        }
        else if (function == 10)
        {
            result = F10(x, xopt, perm, r25, r50, r100, s, w, D);
        }
        else if (function == 11)
        {
            result = F11(x, xopt, perm, r25, r50, r100, s, w, D);
        }
        else if (function == 12)
        {
            result = F12(x, xopt, D);
        }
        else if (function == 13)
        {
            result = F13(x, xopt, perm, r25, r50, r100, s, w, D,overlap);
        }
        //else if (function == 14)
        //{
        //    result = F14(x, xx, perm, r25, r50, r100, s, w, D, overlap);
        //}
        else if (function == 15)
        {
            result = F15(x, xopt,D);
        }

        return result;

    }
}
