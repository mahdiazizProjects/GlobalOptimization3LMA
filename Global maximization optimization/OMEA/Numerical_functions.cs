using System;
using System.Collections.Generic;
using System.Text;
using System.Collections;
namespace OMOA_Namespace
{
    struct Swarms
    {
        public NF[] particles;

    }
    struct Swarms_loc
    {
        public indi[] particles;

    }
    struct BBO
    {
        public double F;
        public double[] x;
        public double Species;
    }
    struct indi
    {
        public double cost;
        public double[] x;
        public double[] y;
        public double[] biasx;
        public double[] biasy;
        public double rho;
        public int visited;
        public bool local_search;
    }
    class Numerical_functions
    {
        public double Normal_Distribution(Random rand, double mean,double stdev)
        {
            //reuse this if you are generating many
            double u1 = rand.NextDouble(); //these are uniform(0,1) random doubles
            double u2 = rand.NextDouble();
            double randStdNormal = Math.Sqrt(-2.0 * Math.Log(u1)) *
                         Math.Sin(2.0 * Math.PI * u2); //random normal(0,1)
          return mean + stdev * randStdNormal;
        }
        public double Cauchy(Random rr, double median,double formFactor)
        {
            double u, v;
            do
            {
                u = 2.0 * rr.NextDouble() - 1.0;
                v = 2.0 * rr.NextDouble() - 1.0;
            }
            while (u * u + v * v > 1.0 || (u == 0.0 && v == 0.0));

            if (u != 0)
                return (median + formFactor * (v / u));
            else
                return (median);
        }

        public int Select_with_Probability(NF[] NFs, double min,int pop_no, Random rr)
        {
            double[] cumsum = cumsum_function(NFs,min, pop_no);
            double rand = rr.NextDouble();
            int result_index = -1;
            for (int ii = 0; ii < pop_no; ii++)
            {
                if (rand > cumsum[ii] && rand <= cumsum[ii + 1])
                {
                    result_index = ii;
                    break;
                }
            }
            if (result_index == -1)
                result_index = pop_no - 1;
            return result_index;
        }
        public int Select_with_Probability(double[] NFs, double min, int pop_no, Random rr)
        {
            double[] cumsum = cumsum_function(NFs, min, pop_no);
            double rand = rr.NextDouble();
            int result_index = -1;
            for (int ii = 0; ii < pop_no; ii++)
            {
                if (rand > cumsum[ii] && rand <= cumsum[ii + 1])
                {
                    result_index = ii;
                    break;
                }
            }
            if (result_index == -1)
                result_index = pop_no - 1;
            return result_index;
        }
        public double[] DE_OP(int index,NF[]POP, Random rr,int N,ArrayList SQ,double F)
        {
            double[] V = new double[N];
            SQ.Remove(index);
            int rand1=rr.Next(SQ.Count);
            int ii1 = Convert.ToInt32(SQ[rand1]);
            SQ.Remove(ii1);
            rand1 = rr.Next(SQ.Count);
            int ii2 = Convert.ToInt32(SQ[rand1]);
            SQ.Remove(ii2);
            for (int ii = 0; ii < N; ii++)
            {
                V[ii] = POP[index].x[ii] + F * (POP[ii1].x[ii] - POP[ii2].x[ii]);
            }
            return V;
        }
        public void PprimtoP(ref NF[] Nodes, ref NF[] NodesP)
        {
            for (int ii = 0; ii < Nodes.Length; ii++)
            {
                Nodes[ii].x =(double[])NodesP[ii].x.Clone();
                Nodes[ii].F = NodesP[ii].F;
            }
        }
        public double[] cumsum_function(NF[] NFs, double min, int pop_no)
        {
            double sum = 0;
            double[] cumsum = new double[pop_no + 1];
            double[] fitness = new double[pop_no];
            for (int ii = 0; ii < pop_no; ii++)
            {
                fitness[ii]=(NFs[ii].F-min);
                sum += fitness[ii];
            }
            for (int ii = 1; ii < pop_no; ii++)
            {
                cumsum[ii] = (fitness[ii] / sum) + cumsum[ii - 1];
            }
            cumsum[pop_no] = 1;
            return cumsum;
        }
        public double[] cumsum_function(double[] NFs, double min, int pop_no)
        {
            double sum = 0;
            double[] cumsum = new double[pop_no + 1];
            double[] fitness = new double[pop_no];
            for (int ii = 0; ii < pop_no; ii++)
            {
                fitness[ii] = (NFs[ii] - min);
                sum += fitness[ii];
            }
            for (int ii = 1; ii < pop_no; ii++)
            {
                cumsum[ii] = (fitness[ii] / sum) + cumsum[ii - 1];
            }
            cumsum[pop_no] = 1;
            return cumsum;
        }
        public NF[] define_pbest(NF[] Nodes, NF[] NeiNodes, int pop_no)
        {
            for (int ii = 0; ii < pop_no; ii++)
            {
                if (NeiNodes[ii].F < Nodes[ii].F)
                {
                    NeiNodes[ii].x = (double[])Nodes[ii].x;
                    NeiNodes[ii].F = Nodes[ii].F;
                }
            }
            return NeiNodes;
        }
        public indi[] define_pbest(indi[] Nodes, indi[] NeiNodes, int pop_no)
        {
            for (int ii = 0; ii < pop_no; ii++)
            {
                if (NeiNodes[ii].cost > Nodes[ii].cost)
                {
                    NeiNodes[ii].x = (double[])Nodes[ii].x;
                    NeiNodes[ii].y = (double[])Nodes[ii].y;
                    NeiNodes[ii].cost = Nodes[ii].cost;
                }
            }
            return NeiNodes;
        }
        public NF Cross_Over_one_point(NF par1, NF par2,int D,int OF,Random rr)
        {
            NF Child = new NF();
            int randp = rr.Next(D);
            double[] temp1 = new double[D];
            double[] temp2 = new double[D];
            for (int ii = 0; ii < randp; ii++)
            {
                temp1[ii] = par1.x[ii];
                temp2[ii] = par2.x[ii];
            }
            for (int ii = randp; ii < D; ii++)
            {
                temp1[ii] = par2.x[ii];
                temp2[ii] = par1.x[ii];
            }
            double f1 = compute_fitness(temp1, OF, D);
            double f2 = compute_fitness(temp2, OF, D);
            if (f1 > f2)
            {
                Child.x = (double[])temp1.Clone();
                Child.F = f1;
            }
            else
            {
                Child.x = (double[])temp2.Clone();
                Child.F = f2;
            }
            return Child;
        }
        public indi Cross_Over_one_point(indi par1, indi par2, int D, Random rr,double[] Anchorx,double[] Anchory,double range,double[] Realx,double[] Realy, double Noise_Factor,double[] gaussian)
        {
            indi Child = new indi();
            int randp = rr.Next(D);
            double[] temp1x = new double[D];
            double[] temp2x = new double[D];
            double[] temp1y = new double[D];
            double[] temp2y = new double[D];
            for (int ii = 0; ii < randp; ii++)
            {
                temp1x[ii] = par1.x[ii];
                temp2x[ii] = par2.x[ii];
                temp1y[ii] = par1.y[ii];
                temp2y[ii] = par2.y[ii];
            }
            for (int ii = randp; ii < D; ii++)
            {
                temp1x[ii] = par2.x[ii];
                temp2x[ii] = par1.x[ii];
                temp1y[ii] = par2.y[ii];
                temp2y[ii] = par1.y[ii];
            }
            double f1 = compute_cost(temp1x, temp1y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);
            double f2 = compute_cost(temp2x, temp2y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);
            if (f1 > f2)
            {
                Child.x = (double[])temp1x.Clone();
                Child.y = (double[])temp1y.Clone();
                Child.cost = f1;
            }
            else
            {
                Child.x = (double[])temp2x.Clone();
                Child.y = (double[])temp2y.Clone();
                Child.cost = f2;
            }
            return Child;
        }
        public NF Cross_Over_two_point(NF par1, NF par2, int D, int OF, Random rr)
        {
            NF Child = new NF();
            int randp1 = rr.Next(D);
            int randp2 = rr.Next(D);
            if (randp1 > randp2)
            {
                int temp = randp1;
                randp1 = randp2;
                randp2 = temp;
            }
            double[] temp1 = new double[D];
            double[] temp2 = new double[D];
            for (int ii = 0; ii < randp1; ii++)
            {
                temp1[ii] = par1.x[ii];
                temp2[ii] = par2.x[ii];
            }
            for (int ii = randp1; ii < randp2; ii++)
            {
                temp1[ii] = par2.x[ii];
                temp2[ii] = par1.x[ii];
            }
            for (int ii = randp2; ii < D; ii++)
            {
                temp1[ii] = par1.x[ii];
                temp2[ii] = par2.x[ii];
            }
            double f1 = compute_fitness(temp1, OF, D);
            double f2 = compute_fitness(temp2, OF, D);
            if (f1 > f2)
            {
                Child.x = (double[])temp1.Clone();
                Child.F = f1;
            }
            else
            {
                Child.x = (double[])temp2.Clone();
                Child.F = f2;
            }
            return Child;
        }
        public NF mask(NF par1, NF par2, int D, int OF, Random rr)
        {
            NF Child = new NF();
            int[] pattern = new int[D];
            for (int ii = 0; ii < D; ii++)
                if (rr.NextDouble() < 0.5)
                    pattern[ii] = 0;
                else
                    pattern[ii] = 1;
            double[] temp1 = new double[D];
            double[] temp2 = new double[D];
            for (int ii = 0; ii < D; ii++)
            {
                if (pattern[ii] == 0)
                {
                    temp1[ii] = par1.x[ii];
                    temp2[ii] = par2.x[ii];
                }
                else 
                {
                    temp1[ii] = par2.x[ii];
                    temp2[ii] = par1.x[ii];
                }
            }
            double f1 = compute_fitness(temp1, OF, D);
            double f2 = compute_fitness(temp2, OF, D);
            if (f1 > f2)
            {
                Child.x = (double[])temp1.Clone();
                Child.F = f1;
            }
            else
            {
                Child.x = (double[])temp2.Clone();
                Child.F = f2;
            }
            return Child;
        }
        public void mutation(ref NF indi, double mu_rate,int D,Random rr)
        {
            for (int ii = 0; ii < D; ii++)
            {
                if (mu_rate > rr.NextDouble())
                {
                    indi.x[ii] = (rr.NextDouble() - 0.5) * 2;
                }
            }
        }
        public void mutation_location(ref indi indi, double mu_rate, int D, Random rr)
        {
            for (int ii = 0; ii < D; ii++)
            {
                if (mu_rate > rr.NextDouble())
                {
                    indi.x[ii] = rr.NextDouble();
                    indi.y[ii] = rr.NextDouble();
                }
            }
        }
        public double[] Mutation_BGA2(ref NF ind, LS[] Pop, int D, Random rr)
        {
            LS indi = new LS();
            indi.x = (double[])ind.x.Clone();
            double num = 16;
            double pai = 1.0 / num;
            int pos = rr.Next(D);
            double rangi = 0; double min = 0; double max = 0;
            int i;
            double dif = 0;
            double sum = 0;
            double[] temp = Min_max_all(pos, Pop);
            min = temp[0];
            max = temp[1];
            rangi = 0.1 * (max - min);

            for (i = 0, dif = 1; i < num; i++, dif /= 2)
            {
                if (rr.NextDouble() < pai)
                {
                    sum += dif;
                }

            }

            double value = indi.x[pos];
            if (sum != 0)
            {
                // Obtain the sign
                if (rr.NextDouble() < 0.5)
                {
                    value += rangi * sum;

                    if (value > max)
                    {
                        value = max;
                    }
                }
                else
                {
                    value -= rangi * sum;

                    if (value < min)
                    {
                        value = min;
                    }
                }
            }

            indi.x[pos] = value;
            return indi.x;
        }
        public double[] Mutation_BGA2(ref double[]x,ref double[]y, indi[] Pop, int D, Random rr)
        {
            indi indi = new indi();
            indi.x = (double[])x.Clone();
            double num = 16;
            double pai = 1.0 / num;
            int pos = rr.Next(D);
            double rangi = 0; double rangi2 = 0; double minx = 0; double maxx = 0; double miny = 0; double maxy = 0;
            int i;
            double dif = 0;
            double sum = 0;
            double[] temp = Min_max_all(pos, Pop);
            minx = temp[0];
            maxx = temp[1];
            rangi = 0.1 * (maxx - minx);
            miny = temp[2];
            maxy = temp[3];
            rangi2 = 0.1 * (maxy - miny);
            for (i = 0, dif = 1; i < num; i++, dif /= 2)
            {
                if (rr.NextDouble() < pai)
                {
                    sum += dif;
                }

            }

            double value = indi.x[pos];
            double value2 = indi.x[pos];
            if (sum != 0)
            {
                // Obtain the sign
                if (rr.NextDouble() < 0.5)
                {
                    value += rangi * sum;

                    if (value > maxx)
                    {
                        value = maxx;
                    }
                    value2 += rangi2 * sum;

                    if (value > maxy)
                    {
                        value = maxy;
                    }
                }
                else
                {
                    value -= rangi * sum;

                    if (value < minx)
                    {
                        value = minx;
                    }
                    value2 -= rangi * sum;

                    if (value < miny)
                    {
                        value = miny;
                    }
                }
            }

            indi.x[pos] = value;
            return indi.x;
        }
        public void Mutation_BGA(ref LS indi, LS[] Pop, int D, Random rr)
        {
            double num = 16;
            double pai = 1.0 / num;
            int pos = rr.Next(D);
            double rangi = 0; double min = 0; double max = 0;
            int i;
            double dif = 0;
            double sum = 0;
            double[] temp = Min_max_all(pos,Pop);
            min = temp[0];
            max = temp[1];
            rangi = 0.1 * (max - min);

            for (i = 0, dif = 1; i < num; i++, dif /= 2)
            {
                if (rr.NextDouble() < pai)
                {
                    sum += dif;
                }

            }

            double value = indi.x[pos];
            if (sum != 0)
            {
                // Obtain the sign
                if (rr.NextDouble() < 0.5)
                {
                    value += rangi * sum;

                    if (value > max)
                    {
                        value = max;
                    }
                }
                else
                {
                    value -= rangi * sum;

                    if (value < min)
                    {
                        value = min;
                    }
                }
            }

            indi.x[pos] = value;
        }
        public double[] Min_max_all(int pos, LS[] Pop)
        {
            double[] temp = new double[2];
            double min = 1000;
            double max = -1000;
            for (int ii = 0; ii < Pop.Length; ii++)
            {
                if (max < Pop[ii].x[pos])
                    max = Pop[ii].x[pos];
                if (min > Pop[ii].x[pos])
                    min = Pop[ii].x[pos];
            }
            temp[0] = min;
            temp[1] = max;
            return temp;
        }
        public double[] Min_max_all(int pos, indi[] Pop)
        {
            double[] temp = new double[4];
            double minx = 1000;
            double maxx = -1000;
            double miny = 1000;
            double maxy= -1000;
            for (int ii = 0; ii < Pop.Length; ii++)
            {
                if (maxx < Pop[ii].x[pos])
                    maxx = Pop[ii].x[pos];
                if (minx > Pop[ii].x[pos])
                    minx = Pop[ii].x[pos];
                if (maxy < Pop[ii].x[pos])
                    maxy = Pop[ii].x[pos];
                if (miny > Pop[ii].x[pos])
                    miny = Pop[ii].x[pos];
            }
            temp[0] = minx;
            temp[1] = maxx;
            temp[2] = miny;
            temp[3] = maxy;
            return temp;
        }
        public double[] Min_max_all(int pos, indi[] Pop,int i)
        {
            double[] temp = new double[2];
            double min = 1000;
            double max = -1000;
            for (int ii = 0; ii < Pop.Length; ii++)
            {
                if (max < Pop[ii].y[pos])
                    max = Pop[ii].y[pos];
                if (min > Pop[ii].y[pos])
                    min = Pop[ii].y[pos];
            }
            temp[0] = min;
            temp[1] = max;
            return temp;
        }
        public double[] generating_jumping_orignal(NF[] Nodes,double[]x ,int pop, int N, int index)
        {
            double min = 1000;
            double max = -100;
            double[] ONodes = new double[N];
            for (int jj = 0; jj < N; jj++)
            {
                for (int ii = 0; ii < pop; ii++)
                {
                    if (min > Nodes[ii].x[jj])
                        min = Nodes[ii].x[jj];
                    if (max < Nodes[ii].x[jj])
                        max = Nodes[ii].x[jj];
                }
            }
            for (int jj = 0; jj < N; jj++)
                ONodes[jj] = min + max - x[jj];
            return ONodes;
        }
        public double[] generating_jumping(LS[] Nodes,double[] Breed ,int pop, int N, int index)
        {
            double min = 1000;
            double max = -100;
            double[] ONodes = new double[N];
            for (int jj = 0; jj < N; jj++)
            {
                for (int ii = 0; ii < pop; ii++)
                {
                    if (min > Nodes[ii].x[jj])
                        min = Nodes[ii].x[jj];
                    if (max < Nodes[ii].x[jj])
                        max = Nodes[ii].x[jj];
                }
            }
            for (int jj = 0; jj < N; jj++)
                ONodes[jj] = min + max - Breed[jj];
            return ONodes;
        }
        public indi generating_jumping(indi[] Nodes, double[] Breed, int pop, int N, int index)
        {
            double minx = 1000;
            double maxx = -100;
            double miny = 1000;
            double maxy = -100;
            double[] ONodesx = new double[N];
            double[] ONodesy = new double[N];
            for (int jj = 0; jj < N; jj++)
            {
                for (int ii = 0; ii < pop; ii++)
                {
                    if (minx> Nodes[ii].x[jj])
                        minx = Nodes[ii].x[jj];
                    if (maxx < Nodes[ii].x[jj])
                        maxx = Nodes[ii].x[jj];
                    if (miny > Nodes[ii].y[jj])
                        miny = Nodes[ii].y[jj];
                    if (maxy < Nodes[ii].y[jj])
                        maxy = Nodes[ii].y[jj];
                }
            }
            for (int jj = 0; jj < N; jj++)
            {
                ONodesx[jj] = minx + maxx - Breed[jj];
                ONodesy[jj] = miny + maxy - Breed[jj];
            }
            indi Result = new indi();
            Result.x = (double[])ONodesx.Clone();
            Result.y = (double[])ONodesy.Clone();
            return Result;
        }
        public int[] sq_rand_gen2(int width, int Max, Random rr)
        {
            int[] sq_rand = new int[width];
            ArrayList squence = new ArrayList();
            for (int ii = 0; ii < Max; ii++)
                squence.Add(ii);
            for (int ii = 0; ii < width; ii++)
            {
                int count = squence.Count;
                int rand = rr.Next(count);
                int _fetch = Convert.ToInt32(squence[rand]);
                squence.Remove(_fetch);
                sq_rand[ii] = _fetch;
            }
            return sq_rand;
        }
        public static void Sort(NF[] X, ref int[] indexes, ref double[] fit)
        {

            for (int ii = 0; ii < X.Length; ii++)
            {
                indexes[ii] = ii;
                fit[ii] = X[ii].F;
            }
            for (int ii = 0; ii < fit.Length - 1; ii++)
                for (int jj = ii + 1; ii < fit.Length; ii++)
                {
                    if (fit[ii] < fit[jj])
                    {
                        double temp = fit[ii];
                        fit[ii] = fit[jj];
                        fit[jj] = temp;
                        int temp2 = indexes[ii];
                        indexes[ii] = indexes[jj];
                        indexes[jj] = temp2;
                    }
                }
        }
        public NF[] mergesort(NF[] CV, int n)
        {
            if (n == 1) return CV;
            int k = n / 2;
            NF[] l2 = new NF[k];
            if ((n % 2) != 0)
                k++;
            NF[] l1 = new NF[k];
            for (int i = 0; i < k; i++)
                l1[i] = CV[i];
            for (int j = k; j < n; j++)
                l2[j - k] = CV[j];
            l1 = mergesort(l1, l1.Length);
            l2 = mergesort(l2, l2.Length);
            return merge(l1, l2);
        }
        private NF[] merge(NF[] a, NF[] b)
        {
            NF[] c = new NF[a.Length + b.Length];
            int k = 0; int p = 0; int q = 0;
            while (a.Length > p && b.Length > q)
            {
                if (a[p].F < b[q].F)
                    c[k++] = b[q++];
                else
                    c[k++] = a[p++];
            }
            while (a.Length > p)
                c[k++] = a[p++];
            while (b.Length > q)
                c[k++] = b[q++];

            return c;
        }
        public indi[] mergesort(indi[] CV, int n)
        {
            if (n == 1) return CV;
            int k = n / 2;
            indi[] l2 = new indi[k];
            if ((n % 2) != 0)
                k++;
            indi[] l1 = new indi[k];
            for (int i = 0; i < k; i++)
                l1[i] = CV[i];
            for (int j = k; j < n; j++)
                l2[j - k] = CV[j];
            l1 = mergesort(l1, l1.Length);
            l2 = mergesort(l2, l2.Length);
            return merge(l1, l2);
        }
        private indi[] merge(indi[] a, indi[] b)
        {
            indi[] c = new indi[a.Length + b.Length];
            int k = 0; int p = 0; int q = 0;
            while (a.Length > p && b.Length > q)
            {
                if (a[p].cost > b[q].cost)
                    c[k++] = b[q++];
                else
                    c[k++] = a[p++];
            }
            while (a.Length > p)
                c[k++] = a[p++];
            while (b.Length > q)
                c[k++] = b[q++];

            return c;
        }
        public BBO[] mergesort(BBO[] CV, int n)
        {
            if (n == 1) return CV;
            int k = n / 2;
            BBO[] l2 = new BBO[k];
            if ((n % 2) != 0)
                k++;
            BBO[] l1 = new BBO[k];
            for (int i = 0; i < k; i++)
                l1[i] = CV[i];
            for (int j = k; j < n; j++)
                l2[j - k] = CV[j];
            l1 = mergesort(l1, l1.Length);
            l2 = mergesort(l2, l2.Length);
            return merge(l1, l2);
        }
        private BBO[] merge(BBO[] a, BBO[] b)
        {
            BBO[] c = new BBO[a.Length + b.Length];
            int k = 0; int p = 0; int q = 0;
            while (a.Length > p && b.Length > q)
            {
                if (a[p].F < b[q].F)
                    c[k++] = b[q++];
                else
                    c[k++] = a[p++];
            }
            while (a.Length > p)
                c[k++] = a[p++];
            while (b.Length > q)
                c[k++] = b[q++];

            return c;
        }
        public LS[] mergesort(LS[] CV, int n)
        {
            if (n == 1) return CV;
            int k = n / 2;
            LS[] l2 = new LS[k];
            if ((n % 2) != 0)
                k++;
            LS[] l1 = new LS[k];
            for (int i = 0; i < k; i++)
                l1[i] = CV[i];
            for (int j = k; j < n; j++)
                l2[j - k] = CV[j];
            l1 = mergesort(l1, l1.Length);
            l2 = mergesort(l2, l2.Length);
            return merge(l1, l2);
        }
        private LS[] merge(LS[] a, LS[] b)
        {
            LS[] c = new LS[a.Length + b.Length];
            int k = 0; int p = 0; int q = 0;
            while (a.Length > p && b.Length > q)
            {
                if (a[p].F < b[q].F)
                    c[k++] = b[q++];
                else
                    c[k++] = a[p++];
            }
            while (a.Length > p)
                c[k++] = a[p++];
            while (b.Length > q)
                c[k++] = b[q++];

            return c;
        }
        public double get_t_quantile(double alpha, int nb_exp, double[,] t_table, bool next = false)
        {
            int pos = -1;
            for (int i = 0; i < 6; i++)
                if (t_table[0, i] == alpha)
                {
                    pos = i;
                    break;
                }
            if (next) pos++; //(1-alpha/2) quantile
            if (nb_exp > 100) nb_exp = 101;
            return (t_table[nb_exp, pos]);
        }
        public double compute_fitness(double[] X, int objective_function, int S)
        {
            double result = 0;
            double[] x = (double[])X.Clone();
            #region functions
            if (objective_function == 1)
            {
                for (int ii = 0; ii < S; ii++)
                {
                    x[ii] = x[ii] * 500;
                    result += x[ii] * Math.Sin(Math.Sqrt(Math.Abs(x[ii])));
                }
            }
            else if (objective_function == 2)
            {
                for (int ii = 0; ii < S; ii++)
                {
                    x[ii] = x[ii] * 5.12;
                    result += Math.Pow(x[ii], 2) - 10 * Math.Cos(2 * Math.PI * x[ii]) + 10;
                }
                result = -result;
            }
            else if (objective_function == 3)
            {

                double mean1 = 0;
                double mean2 = 0;
                for (int ii = 0; ii < S; ii++)
                {
                    x[ii] = x[ii] * 32;
                    mean1 += Math.Cos(2 * Math.PI * x[ii]);
                    mean2 += Math.Pow(x[ii], 2);
                }
                mean1 = mean1 / S;
                mean2 = mean2 / S;
                result = (20 * (Math.Exp(-0.2 * Math.Sqrt(mean2))) + Math.Exp(mean1) - 20-Math.E);
            }
            else if (objective_function == 4)
            {
                double sum = 0;
                double prod = 1;
                for (int ii = 1; ii < S; ii++)
                {
                    x[ii-1] = x[ii-1] * 600;
                    sum += Math.Pow(x[ii-1], 2);
                    prod *= Math.Cos(x[ii-1] / Math.Sqrt(Convert.ToDouble(ii)));
                }
                result =1- (sum / 4000) - prod;
            }
            else if (objective_function == 5)
            {
                double[] y = new double[S];
                double u=0; 
                double res_y=0;
                for (int ii = 0; ii < S; ii++)
                {
                    x[ii] = 50 * x[ii];
                    y[ii] = 1 + (x[ii] + 1) / 4;
                    u+=compute_u(x[ii],10,100,1);
                   
                }
                for (int ii = 0; ii < S-1; ii++)
                {
                    res_y += Math.Pow((y[ii] - 1), 2) * (1 + 10 * Math.Pow(Math.Sin(Math.PI * y[ii + 1]), 2));
                }
                result=-(Math.PI / S)*(10 * Math.Pow(Math.Sin(Math.PI * y[0]),2) - res_y-Math.Pow(y[S-1],2)) - u;
                
            }
            else if (objective_function == 6)
            {
                double res_y = 0;
                double u = 0;
                for (int ii = 0; ii < S; ii++)
                {
                    x[ii] = 50 * x[ii];
                    u += compute_u(x[ii], 5, 100, 1);
                }
                for (int ii = 0; ii < S-1; ii++)
                {
                    res_y += Math.Pow((x[ii] - 1), 2) * (1 + Math.Pow(Math.Sin(3 * Math.PI * x[ii + 1]), 2));
                }
                 result = -(0.1)*( Math.Pow(Math.Sin(3*Math.PI * x[0]), 2) + res_y + Math.Pow(x[S - 1], 2)*(1+Math.Pow(Math.Sin(2*Math.PI*x[S-1]),2)))-u;
            }
            else if (objective_function == 7)
            {
                for (int ii = 0; ii < S; ii++)
                    x[ii] = Math.PI * x[ii];
                for (int ii = 0; ii < S; ii++)
                    result += Math.Sin(x[ii]) * Math.Pow(Math.Sin((ii * Math.Pow(x[ii], 2)) / Math.PI), 2 * S);
            }
            else if (objective_function == 8)
            {
                double h = 0;
                double g = 0;
                for (int ii = 0; ii < S; ii++)
                {
                    g = Math.Pow(Math.Sin(5 * Math.PI * x[ii] + 0.5), 2);
                    h = Math.Exp(-2 * Math.Log10(2) * (Math.Pow(x[ii] - 0.1, 2) / 0.64));
                    result += g * h;
                }
            }
            else if (objective_function == 9)
            {
                for (int ii = 0; ii < S; ii++)
                {
                    x[ii] = 100 * x[ii];
                    result += Math.Pow(x[ii], 2);
                }
                result = -result;
            }
            else if (objective_function == 10)
            {
                double x1 = 0;
                double x2 = 1;
                for (int ii = 0; ii < S; ii++)
                {
                    x[ii] = 10 * x[ii];
                    x1 += Math.Abs(x[ii]);
                    x2 *= Math.Abs(x[ii]);
                }
                result = -x1 - x2;
            }
            else if (objective_function == 11)
            {
                double max = -1000000;
                for (int ii = 0; ii < S; ii++)
                {
                    x[ii] = 100 * x[ii];
                    max = Math.Max(Math.Abs(x[ii]), max);
                }
                result = -max;
            }
            else if (objective_function == 12)
            {
                double rand = 0;
                for (int ii = 0; ii < S; ii++)
                {
                    x[ii] = 10 * x[ii];
                    result += (ii) * Math.Pow(x[ii], 4);
                }
                Random rr = new Random();
                rand = Normal_Distribution(rr, 0, 1);
                result = -result - rand;
            }
            else if (objective_function == 13)
            {
                for (int ii = 0; ii < S; ii++)
                {
                    x[ii] = 2 * x[ii];
                }
                for (int ii = 0; ii < S - 1; ii++)
                {
                    result += 100 * (Math.Pow((x[ii + 1] - Math.Pow(x[ii], 2)), 2)) + Math.Pow(1 - x[ii], 2);
                }
                result = -result;
            }

            else if (objective_function == 14)
            {
                for (int ii = 0; ii < S; ii++)
                {
                    x[ii] = 4 * x[ii];
                }
                double min=100000000;
                for (double jj = 1; jj < S; jj++)
                {
                    double res = 0;
                    for (int ii = 0; ii < S; ii++)
                    {

                        res += Math.Pow((1 / (1 + Math.Exp(-x[ii]))), 2);
                    }
                    res += (Math.Pow(jj - 1, 0.15)) / 15;
                    min = Math.Min(res, min);
                }
                result = min; 
            }
            else if (objective_function == 15)
            {
                for (int ii = 0; ii < S; ii++)
                {
                    x[ii] = 10 * x[ii];
                        result += Math.Pow(1000000, ((ii) / S - 1)) * Math.Pow(x[ii], 2);
                }
                result = -result;

            }
            else if (objective_function == 22)
            {
                double temp = 0;
                for (int ii = 0; ii < S; ii++)
                {
                    x[ii] = 65.536 * x[ii];
                }
                for (int ii = 0; ii < S; ii++)
                {
                    for (int jj = 0; jj < ii; jj++)
                    {
                        temp += x[ii];
                    }
                    result += Math.Pow(temp, 2);
                }
                result = -result;
            }
#endregion
            return result;
        }

        public double compute_fitness(double[] X, int objective_function, int S,int[]Perm,double[,]Orthogonal_Matrix,int M,double[]o )
        {
           
            double result = 0;
            double[] x = (double[])X.Clone();
            #region functions

            if (objective_function == 16)
            {
                double[] z = new double[S];
                for (int ii = 0; ii < S; ii++)
                {
                    x[ii] = 10 * x[ii];
                    double temp = 0;
                    for (int jj = 0; jj < M; jj++)
                        temp += x[jj] * Orthogonal_Matrix[ii, jj]*10;
                    z[ii] = temp;
                    result += Math.Pow(1000000, ((ii) / S - 1)) * Math.Pow(z[ii], 2);
                }
                result = -result;
            }
            else if (objective_function == 17)
            {
                double[] z = new double[M];
                for (int ii = 0; ii < S; ii++)
                {
                    x[ii] = x[ii] * 10;
                    x[ii] = x[ii] - o[ii] * 10;
                }
                for (int ii = 0; ii < M; ii++)
                {
                    double temp = 0;
                    for (int jj = 0; jj < M; jj++)
                        temp += x[Perm[jj]] * Orthogonal_Matrix[ii, jj]*10;
                    z[ii] = temp;
                        result += Math.Pow(1000000, ((ii ) / S - 1)) * Math.Pow(z[ii], 2);
                }

                for (int ii = M; ii < S; ii++)
                {
                        result += Math.Pow(1000000, ((ii) / S - 1)) * Math.Pow(x[Perm[ii]], 2);  
                }
                result = -result;
            }
            else if (objective_function == 18)
            {

                double[] z = new double[M];
                for (int ii = 0; ii < S; ii++)
                {
                    x[ii] = x[ii] * 5.12;
                    x[ii] = x[ii] - o[ii]*5.12;
                }
                for (int ii = 0; ii < M; ii++)
                {
                    double temp = 0;
                    for (int jj = 0; jj < M; jj++)
                        temp += x[Perm[jj]] * Orthogonal_Matrix[ii, jj]*5.12;
                    z[ii] = temp;
                    result += Math.Pow(z[ii], 2) - (10 * Math.Cos(2 * Math.PI * z[ii])) + 10;
                }
                for (int ii = M; ii < S; ii++)
                    result += Math.Pow(x[Perm[ii]], 2) - (10 * Math.Cos(2 * Math.PI * x[Perm[ii]])) + 10;
                result = -result;
            }
            else if (objective_function == 19)
            {
                double[] z = new double[M];
                double mean1 = 0;
                double mean2 = 0;
                for (int ii = 0; ii < S; ii++)
                {
                    x[ii] = x[ii] * 32;
                    x[ii] = x[ii] - o[ii] * 32;
                }
                for (int ii = 0; ii < M; ii++)
                {
                    double temp = 0;
                    for (int jj = 0; jj < M; jj++)
                        temp += x[Perm[jj]] * Orthogonal_Matrix[ii, jj]*32;
                    z[ii] = temp;
                    mean1 += Math.Cos(2 * Math.PI * z[ii]);
                    mean2 += Math.Pow(z[ii], 2);
                }
                mean1 = mean1 / M;
                mean2 = mean2 / M;
                result = -(20 * (Math.Exp(-0.2 * Math.Sqrt(mean2))) + Math.Exp(mean1) - 20 - Math.E);


                mean1 = 0;
                mean2 = 0;
                for (int ii = M; ii < S; ii++)
                {
                    mean1 += Math.Cos(2 * Math.PI * x[Perm[ii]]);
                    mean2 += Math.Pow(x[Perm[ii]], 2);
                }
                mean1 = mean1 / (S - M);
                mean2 = mean2 / (S - M);
                result += -(20 * (Math.Exp(-0.2 * Math.Sqrt(mean2))) + Math.Exp(mean1) - 20 - Math.E);
            }
            else if (objective_function == 20)
            {
                double max = -1000000;
                double[] z = new double[M];
                for (int ii = 0; ii < S; ii++)
                {
                    x[ii] = 100 * x[ii];
                    x[ii] = x[ii] - o[ii] * 100;
                }

                for (int ii = 0; ii < M; ii++)
                {
                    double temp = 0;
                    for (int jj = 0; jj < M; jj++)
                        temp += x[Perm[jj]] * Orthogonal_Matrix[ii, jj]*100;
                    z[ii] = temp;
                    max = Math.Max(Math.Abs(z[ii]), max);
                }
                result = max;

                max = -1000000;
                for (int ii = M; ii < S; ii++)
                {
                    max = Math.Max(Math.Abs(x[Perm[ii]]), max);
                }
                result =-result- max;
            }
            else if (objective_function == 21)
            {

                double[] z = new double[M];
                for (int ii = 0; ii < S; ii++)
                {
                    x[ii] = 2 * x[ii];
                    x[ii] = x[ii] - o[ii] * 2;
                }
                for (int ii = 0; ii < M - 1; ii++)
                {
                    double temp = 0;
                    for (int jj = 0; jj < M; jj++)
                        temp += x[Perm[jj]] * Orthogonal_Matrix[ii, jj]*2;
                    z[ii] = temp;
                    result += 100 * (Math.Pow((z[ii + 1] - Math.Pow(z[ii], 2)), 2)) + Math.Pow(1 - z[ii], 2);
                }

                for (int ii = M - 1; ii < S - 1; ii++)
                {
                    result += 100 * (Math.Pow((x[Perm[ii+1]] - Math.Pow(x[Perm[ii]], 2)), 2)) + Math.Pow(1 - x[Perm[ii]], 2);
                }
                result = -result;
            }
            else if(objective_function==22)
            {
                 //F10: D/2m-group Shifted and m-rotated Rastrigin’s Function
            
                int G = S / M / 2;
                double[] Y = new double[M];
                for (int ii = 0; ii < S; ii++)
                {
                    x[ii] = x[ii] - o[ii];
                    x[ii] = 5 * x[ii];
                }
                for (int ii = 0; ii < G; ii++)
                {
                    for (int jj = (ii * M); jj < ((ii+1) * M); jj++)
                        Y[jj] = x[Perm[jj]];
                    result += The_Rotated_Rastrigins_Function(Orthogonal_Matrix, Y, M);
                }
                double[] Z = new double[M * G];
                for (int ii = 0; ii < Z.Length; ii++)
                    Z[ii] = x[Perm[ii]];
                result += Rastrigin_Function(Z, M * G);
            
            }
            //F11: D/2m-group Shifted and m-rotated Ackley’s Function
            else if (objective_function == 23)
            {
                int G = S / M / 2;
                double[] Y = new double[M];
                for (int ii = 0; ii < S; ii++)
                {
                    x[ii] = x[ii] - o[ii];
                    x[ii] = 32 * x[ii];
                }
                for (int ii = 0; ii < G; ii++)
                {
                    int k = 0;
                    for (int jj = (ii * M); jj < ((ii + 1) * M); jj++)
                        Y[k++] = x[Perm[jj]];
                     
                    result += The_Rotated_Ackley_Function(Orthogonal_Matrix, Y, M);
                }
                double[] Z = new double[M * G];
                for (int ii = 0; ii < Z.Length; ii++)
                    Z[ii] = x[Perm[ii]];
                result += Ackley_Function(Z, M * G);
            }
            ////F13: D/2m-group Shifted m-dimensional Rosenbrock’s Function
            else if (objective_function == 24)
            {
                int G = S / M / 2;
                double[] Y = new double[M];
                for (int ii = 0; ii < S; ii++)
                {

                    x[ii] = x[ii] - o[ii];
                    x[ii] = 100 * x[ii];
                }
                for (int ii = 0; ii < G; ii++)
                {
                    int kk = 0;
                    for (int jj = (ii * M ); jj < ((ii+1) * M); jj++)
                        Y[kk++] = x[Perm[jj]];
                    result += Rosenbrock_Function(Y, M);
                }
                double[] Z = new double[M * G];
                for (int ii = 0; ii < Z.Length; ii++)
                    Z[ii] = x[Perm[ii]];
                result += Sphere_Function(Z, M * G);
            }
            ////F15: D/m-group Shifted and m-rotated Rastrigin’s Function
            else if (objective_function == 25)
            {
                int G = S / M;
                double[] Y = new double[M];
                for (int ii = 0; ii < S; ii++)
                {

                    x[ii] = x[ii] - o[ii];
                    x[ii] = 5 * x[ii];
                }
                for (int ii = 0; ii < G; ii++)
                {
                    int k = 0;
                    for (int jj = ((ii) * M ); jj < ((ii+1) * M); jj++)
                        Y[k++] = x[jj];
                    result += The_Rotated_Rastrigins_Function(Orthogonal_Matrix, Y, M);
                }
            }
            //////////F16: D/m-group Shifted and m-rotated Ackley’s Function
            else if (objective_function == 26)
            {
                int G = S / M;
                double[] Y = new double[M];
                for (int ii = 0; ii < S; ii++)
                {
                    x[ii] = x[ii] - o[ii];
                    x[ii] = 32 * x[ii];
                }
                for (int ii = 0; ii < G; ii++)
                {
                    int k = 0;
                    for (int jj = (ii * M); jj < ((ii+1) * M); jj++)
                        Y[k++] = x[jj];
                    result += The_Rotated_Ackley_Function(Orthogonal_Matrix, Y, M);
                }
            }
            //////F18: D/m-group Shifted m-dimensional Rosenbrock’s Function
            else if (objective_function == 27)
            {
                int G = S / M;
                double[] Y = new double[M];
                for (int ii = 0; ii < S; ii++)
                {

                    x[ii] = x[ii] - o[ii];
                    x[ii] = 100 * x[ii];
                }
                for (int ii = 0; ii < G; ii++)
                {
                    int k = 0;
                    for (int jj = ((ii ) * M ); jj < ((ii+1) * M); jj++)
                        Y[k++] = x[jj];
                    result += Rosenbrock_Function(Y, M);
                }
            }
            
            #endregion
            return result;
        }
        public double Ackley_Function(double[] x, int S)
        {
            double result = 0;
            double mean1 = 0;
            double mean2 = 0;
            for (int ii = 0; ii < S; ii++)
            {
                mean1 += Math.Cos(2 * Math.PI * x[ii]);
                mean2 += Math.Pow(x[ii], 2);
            }
            mean1 = mean1 / S;
            mean2 = mean2 / S;
            result = 20 * (Math.Exp(-0.2 * Math.Sqrt(mean2))) + Math.Exp(mean1) - 20 - Math.E;
            return result;
        }
        public double Rastrigin_Function(double[] x, int S)
        {
            double result = 0;
            for (int ii = 0; ii < S; ii++)
            {
                result += Math.Pow(x[ii], 2) - 10 * Math.Cos(2 * Math.PI * x[ii]) + 10;
            }
            return -result;
        }
        public double Sphere_Function(double[] x, int S)
        {
            double result = 0;
            for (int ii = 0; ii < S; ii++)
            {
                x[ii] = 100 * x[ii];
                result += Math.Pow(x[ii], 2);
            }
            return -result;
        }
        public double Rosenbrock_Function(double[] x, int S)
        {
            double result = 0;
            for (int ii = 0; ii < S - 1; ii++)
            {
                result += 100 * (Math.Pow((Math.Pow(x[ii], 2) - x[ii + 1]), 2)) - Math.Pow(1-x[ii], 2);
            }
            return -result;
        }
        public double The_Rotated_Ackley_Function(double[,] Orthogonal_Matrix, double[] x, int S)
        {
            double mean1 = 0; double mean2 = 0; double result = 0; double temp = 0; double[] z = new double[S];
            for (int ii = 0; ii < S; ii++)
            {
                temp = 0;
                for (int jj = 0; jj < S; jj++)
                    temp += x[jj] * Orthogonal_Matrix[jj, ii];
                z[ii] = temp;
                mean1 += Math.Cos(2 * Math.PI * z[ii]);
                mean2 += Math.Pow(z[ii], 2);
            }
            mean1 = mean1 / S;
            mean2 = mean2 / S;
            result = 20 * (Math.Exp(-0.2 * Math.Sqrt(mean2))) + Math.Exp(mean1) - 20 - Math.E;
            return result;
        }
        public double The_Rotated_Rastrigins_Function(double[,] Orthogonal_Matrix, double[] x, int S)
        {
            double result = 0;
            double[] z = new double[S];
            for (int ii = 0; ii < S; ii++)
            {
                double temp = 0;
                for (int jj = 0; jj < S; jj++)
                    temp += x[jj] * Orthogonal_Matrix[jj, ii];
                z[ii] = temp;
                result += Math.Pow(z[ii], 2) - 10 * Math.Cos(2 * Math.PI * z[ii]) + 10;
            }
            return -result;
        }
        public double[] node_constructing(int N, Random rr)
        {
            double[] x = new double[N];
            for (int ii = 0; ii < N; ii++)
            {
                int temp = rr.Next(2);
                double temp2=rr.NextDouble();
                if (temp == 1)
                    temp2 = -temp2;
                x[ii] = temp2;
            }
            return x;
        }
        public double[] node_constructing2(int N, Random rr,int u,int l)
        {
            double[] x = new double[N];
            for (int ii = 0; ii < N; ii++)
            {
                x[ii] = l + 2 * u * rr.NextDouble();
            }
            return x;
        }
        public double[] Initialization_Hypercube(double[] Elitex, int N,Random rr, double delta)
        {
            double[] x = new double[N];
            for (int ii = 0; ii < N; ii++)
            {
                int temp = rr.Next(2);
                double temp2 = rr.NextDouble();
                if (temp == 1)
                    temp2 = -temp2;
                x[ii] = Elitex[ii]+ temp2* delta;
            }
            return x;
        }
        public double[] Tau_constructing_ACO(int N, Random rr,double tau0)
        {
            double[] tau = new double[N];
            for (int ii = 0; ii < N; ii++)
            {
                tau[ii] = tau0 * (1);
            }
            return tau;
        }
        public void Tau_update(ref double[] tau, double Gdr)
        {
            for (int ii = 0; ii < tau.Length; ii++)
                tau[ii] = (1 - Gdr) * tau[ii];
        }
        public double compute_u(double x, double a, double k, double p)
        {
            if (x > a)
                return k * Math.Pow((x - a), p);
            else if (x < -a)
                return k * Math.Pow((-x - a), p);
            else
                return 0;
        }
        public int[] Structures(int ii, int pop_no, string structure)
        {
            int[] neighbours = { 0, 1 };
            if (structure == "Ring")
                neighbours = Ring(ii, pop_no);
            else if (structure == "Star")
                neighbours = Star_structure(ii, pop_no);
            else if (structure == "Btree")
                neighbours = Btree_structure(ii, pop_no);
            else if (structure == "Bipartie")
                neighbours = Bipartie_structure(ii, pop_no);
            else if (structure == "Cluster")
                neighbours = Cluster_structure(ii, pop_no, 5);
            else if (structure == "Grid")
                neighbours = grid_structure(ii, pop_no);
            else if (structure == "Ladder")
                neighbours = ladder_structure(ii, pop_no);
            else if (structure == "Cross-Ladder")
                neighbours = cross_ladder_structure(ii, pop_no);
            else if (structure == "Cellular")
                neighbours = cellular_structure(ii, pop_no);
            return neighbours;
        }
        public int[] cellular_structure(int ii, int pop_no)
        {
            double y = 0.5;
            int L = Convert.ToInt32(Math.Floor(Math.Pow(Convert.ToDouble(pop_no), y)));
            int row = 0; int column = 0;
            int pop_no2 = Convert.ToInt32(Math.Pow(L, 2));
            int[] neighbor = new int[4];
            row = Convert.ToInt32(Math.Ceiling(Convert.ToDouble(ii / L)));
            column = ((ii) % L);
            neighbor[0] = ii - 1;
            neighbor[1] = ii + 1;
            neighbor[2] = ii + L;
            neighbor[3] = ii - L;
            if (column == 0)
                neighbor[0] = ii + L - 1;
            if (column == L - 1)
                neighbor[1] = ii - L + 1;
            if (row == L - 1 && pop_no == pop_no2)
                neighbor[2] = ii - (L - 1) * L;
            else if (row == L && pop_no2 != pop_no)
                neighbor[2] = ii - (L) * L;
            if (row == 0 && pop_no == pop_no2)
                neighbor[3] = ii + (L - 1) * L;
            else if (row == 0 && pop_no2 != pop_no)
                neighbor[3] = ii + (L) * L;
            if (neighbor[3] >= pop_no)
                neighbor[3] = pop_no-1;
            if (neighbor[2] >= pop_no)
                neighbor[2] = pop_no - 1;
            if (neighbor[1] >= pop_no)
                neighbor[1] = pop_no - 1;
            if (neighbor[0] >= pop_no)
                neighbor[0] = pop_no - 1;

            return neighbor;
        }
        public int[] grid_structure(int ii, int pop_no)
        {
            double y = 0.5;
            int L = Convert.ToInt32(Math.Floor(Math.Pow(Convert.ToDouble(pop_no), y)));
            int row = 0; int column = 0;
            int pop_no2 = Convert.ToInt32(Math.Pow(L, 2));
            ArrayList neighbour = new ArrayList();
            row = Convert.ToInt32(Math.Ceiling(Convert.ToDouble(ii / L)));
            column = ((ii) % L);
            neighbour.Add(ii - 1);
            neighbour.Add(ii + 1);
            neighbour.Add(ii + L);
            neighbour.Add(ii - L);
            //if (column == 0)
            //    neighbour.Remove(ii - 1);
            //if (column == L - 1)
            //    neighbour.Remove(ii + 1);
            //if (row == L - 1 && pop_no == pop_no2)
            //    neighbour.Remove(ii + L);
            //else if (row == L && pop_no2 != pop_no)
            //    neighbour.Remove(ii + L);
            //if (row == 0 && pop_no == pop_no2)
            //    neighbour.Remove(ii - L);
            //else if (row == 0 && pop_no2 != pop_no)
            //    neighbour.Remove(ii - L);
            for (int jj = 0; jj < neighbour.Count; jj++)
            {
                int temp = Convert.ToInt32(neighbour[jj]);
                if (temp >= pop_no || temp < 0)
                {
                    neighbour.RemoveAt(jj);
                    jj--;
                }
            }

            int[] neighbor = (int[])neighbour.ToArray(Type.GetType("System.Int32"));
            return neighbor;
        }
        public int[] ladder_structure(int ii, int pop_no)
        {
            double pop = Convert.ToDouble(pop_no);
            int L = Convert.ToInt32(Math.Floor(pop / 2));
            ArrayList neighbour = new ArrayList();
            neighbour.Add(ii - 1);
            neighbour.Add(ii + 1);
            if (ii < L)
            {
                if (neighbour.Contains(-1))
                    neighbour.Remove(ii - 1);
                if (neighbour.Contains(L))
                    neighbour.Remove(L);
                neighbour.Add(ii + L);
            }
            else
            {
                if (neighbour.Contains(L))
                    neighbour.Remove(L);
                if (neighbour.Contains(pop_no))
                    neighbour.Remove(pop_no);
                neighbour.Add(ii - L);
            }
            int[] neighbor = (int[])neighbour.ToArray(Type.GetType("System.Int32"));
            return neighbor;
        }
        public int[] cross_ladder_structure(int ii, int pop_no)
        {
            double pop = Convert.ToDouble(pop_no);
            int L = Convert.ToInt32(Math.Floor(pop / 2));
            ArrayList neighbour = new ArrayList();
            neighbour.Add(ii - 1);
            neighbour.Add(ii + 1);
            if (ii < L)
            {
                if (neighbour.Contains(-1))
                    neighbour.Remove(-1);
                if (neighbour.Contains(L))
                    neighbour.Remove(L);
                neighbour.Add(ii + L);
                neighbour.Add(ii + L - 1);
                if (neighbour.Contains(L))
                    neighbour.Remove(L);
                neighbour.Add(ii + L + 1);
                if (neighbour.Contains(pop_no))
                    neighbour.Remove(pop_no);
            }
            else
            {
                if (neighbour.Contains(L))
                    neighbour.Remove(L);
                if (neighbour.Contains(pop_no))
                    neighbour.Remove(pop_no);
                neighbour.Add(ii - L);
                neighbour.Add(ii - L - 1);
                if (neighbour.Contains(-1))
                    neighbour.Remove(-1);
                neighbour.Add(ii - L + 1);
                if (neighbour.Contains(L))
                    neighbour.Remove(L);
            }
            int[] neighbor = (int[])neighbour.ToArray(Type.GetType("System.Int32"));
            return neighbor;
        }
        public int[] Btree_structure(int ii, int pop_no)
        {
            ArrayList neighbour = new ArrayList();
            Queue<int> Q = new Queue<int>();
            do
            {
                if ((2 * ii + 1) < pop_no)
                {
                    neighbour.Add(2 * ii + 1);
                    Q.Enqueue(2 * ii + 1);
                }
                if ((2 * ii + 2) < pop_no)
                {
                    neighbour.Add(2 * ii + 2);
                    Q.Enqueue(2 * ii + 2);
                }
                if (Q.Count > 0)
                {
                    ii = Q.Dequeue();
                }
            } while (2 * ii + 1 < pop_no || Q.Count > 0);
            int[] neighbor = (int[])neighbour.ToArray(Type.GetType("System.Int32"));
            return neighbor;
        }
        public int[] Star_structure(int ii, int pop_no)
        {
            int[] neighbor = new int[pop_no - 1];
            int kk = 0;
            for (int jj = 0; jj < pop_no; jj++)
                if (jj != ii)
                    neighbor[kk++] = jj;

            return neighbor;
        }
        public int[] Cluster_structure(int ii, int pop_no, int cluster_size)
        {
            ArrayList neighbors = new ArrayList();
            int first_neighbour = (ii / cluster_size) * cluster_size;
            int[] neighbour1 = new int[2];
            int end_neighbour = Convert.ToInt32(Math.Ceiling(Convert.ToDouble(ii + 1) / cluster_size)) * cluster_size;
            if ((ii % cluster_size) != 0 || ii == 0)
            {
                for (int jj = first_neighbour; jj < end_neighbour; jj++)
                    neighbors.Add(jj);
            }
            else if ((ii % 5) == 0 && ii != 0)
            {
                int cluster = ii / (cluster_size);
                neighbour1[0] = (cluster - 1) * 5;
                if (neighbour1[0] <= 0)
                    neighbour1[0] = cluster_size;
                neighbour1[1] = (cluster + 1) * 5;
                if (neighbour1[1] >= pop_no)
                    neighbour1[1] = pop_no - 1;
                neighbors.Add(neighbour1[0]); neighbors.Add(neighbour1[1]);
                first_neighbour = ((ii - 1) / cluster_size) * cluster_size + 1;
                end_neighbour = Convert.ToInt32(Math.Ceiling(Convert.ToDouble(ii - 1) / cluster_size)) * cluster_size;
                for (int jj = first_neighbour; jj < end_neighbour; jj++)
                    if (!neighbors.Contains(jj))
                        neighbors.Add(jj);
            }
            int[] neighbor = (int[])neighbors.ToArray(System.Type.GetType("System.Int32"));

            return neighbor;
        }
        public int[] Bipartie_structure(int ii, int pop_no)
        {
            double pop = Convert.ToDouble(pop_no);
            int L = Convert.ToInt32(Math.Floor(pop / 2));
            int[] neighbor = new int[L];
            if (ii < L)
            {
                for (int jj = L, kk = 0; jj < pop_no; jj++, kk++)
                    neighbor[kk] = jj;
            }
            else
            {
                for (int jj = 0; jj < L; jj++)
                    neighbor[jj] = jj;
            }
            return neighbor;
        }
        public int[] Ring(int ii, int pop_no)
        {
            int[] neighbor = new int[2];
            neighbor[0] = ii - 1;
            neighbor[1] = ii + 1;
            if (neighbor[0] < 0)
                neighbor[0] = pop_no - 1;
            if (neighbor[1] > (pop_no - 1))
                neighbor[1] = 0;
            return neighbor;
        }
        public int[] Ring_withmyself(int ii, int pop_no)
        {
            int[] neighbor = new int[3];
            neighbor[0] = ii - 1;
            neighbor[1] = ii + 1;
            neighbor[2] = ii;
            if (neighbor[0] < 0)
                neighbor[0] = pop_no - 1;
            if (neighbor[1] > (pop_no - 1))
                neighbor[1] = 0;
            return neighbor;
        }
        public double[] Opposition_Learning(double[] x)
        {
            double[]Ox=new double[x.Length];
            for (int ii = 0; ii < x.Length; ii++)
                Ox[ii] = -x[ii];
            return Ox;
        }
        public void SolisWets(ref double[] Sol,ref double fitness_Sol, int N, ref  double[] bias, ref double rho, int MaxEval, Random rand, int OF,int[]Perm,double[,]Orthogonal_Matrix,int M,double[] o)
        {
            int numSuccess = 0;
            int NumFailed = 0;
            double[]dif=new double[N];
            double[]Solprim=new double[N];
            double[] Solzegon = new double[N];
            double fitness_Sol_prim = 0;
            double fitness_Sol_zegon = 0;
            for (int numEval = 0; numEval < MaxEval; numEval++)
            {
                for (int ii = 0; ii < N; ii++)
                {
                    dif[ii] = Normal_Distribution(rand, 0, rho);
                    Solprim[ii] = Sol[ii] + bias[ii] + dif[ii];
                }
                Solprim = Node_treatment(Solprim, N, rand);
                if (OF < 16 || OF == 22)
                    fitness_Sol_prim = compute_fitness(Solprim, OF, N);
                else
                    fitness_Sol_prim = compute_fitness(Solprim, OF, N, Perm, Orthogonal_Matrix, M, o);
                if (fitness_Sol < fitness_Sol_prim)
                {
                    Sol = (double[])Solprim.Clone();
                    fitness_Sol = fitness_Sol_prim;
                    for (int ii = 0; ii < N; ii++)
                    {
                        bias[ii] = 0.2 * bias[ii] + 0.4 * (dif[ii] + bias[ii]);
                    }
                    numSuccess++;
                    NumFailed = 0;
                }
                else
                {
                    for (int ii = 0; ii < N; ii++)
                        Solzegon[ii] = Sol[ii] - bias[ii] - dif[ii];
                    Solzegon=Node_treatment(Solzegon, N, rand);
                    if (OF < 16 || OF == 22)
                        fitness_Sol_zegon = compute_fitness(Solzegon, OF, N);
                    else
                        fitness_Sol_zegon = compute_fitness(Solzegon, OF, N, Perm, Orthogonal_Matrix, M, o);
                    if (fitness_Sol_zegon > fitness_Sol)
                    {
                        Sol = (double[])Solzegon.Clone();
                        fitness_Sol = fitness_Sol_zegon;
                        for (int ii = 0; ii < N; ii++)
                        {
                            bias[ii] = bias[ii] - 0.4 * (dif[ii] + bias[ii]);
                        }
                        numSuccess++;
                        NumFailed = 0;

                    }
                    else
                    {
                        NumFailed++;
                        numSuccess = 0;
                    }

                }

                if (numSuccess > 5)
                {
                    rho *= 2;
                    numSuccess = 0;
                }
                else if (NumFailed > 3)
                {
                    rho /= 2;
                    NumFailed = 0;
                }
            }

        }
        public void SolisWets(ref double[] Solx, ref double[] Soly, ref double fitness_Sol, int N, ref  double[] biasx,ref double[] biasy, ref double rho, int MaxEval, Random rand, double range, double[] rx, double[] ry, double[] gaussian, double noise_factor, double[] Anchorx, double[] Anchory)
        {
            int numSuccess = 0;
            int NumFailed = 0;
            double[] difx = new double[N];
            double[] dify = new double[N];
            double[] Solprimx = new double[N];
            double[] Solprimy = new double[N];
            double[] Solzegonx = new double[N];
            double[] Solzegony = new double[N];
            double fitness_Sol_prim = 0;
            double fitness_Sol_zegon = 0;
            for (int numEval = 0; numEval < MaxEval; numEval++)
            {
                for (int ii = 0; ii < N; ii++)
                {
                    difx[ii] = Normal_Distribution(rand, 0, rho);
                    Solprimx[ii] = Solx[ii] + biasx[ii] + difx[ii];
                    dify[ii] = Normal_Distribution(rand, 0, rho);
                    Solprimy[ii] = Soly[ii] + biasy[ii] + dify[ii];
                }
                Node_treatment_localization(ref Solprimx, ref Solprimy, N, rand);
                fitness_Sol_prim = compute_cost(Solprimx, Solprimy, Anchorx, Anchory, range, range, rx, ry, noise_factor, gaussian);
                if (fitness_Sol > fitness_Sol_prim)
                {
                    Solx = (double[])Solprimx.Clone();
                    Soly = (double[])Solprimy.Clone();
                    fitness_Sol = fitness_Sol_prim;
                    for (int ii = 0; ii < N; ii++)
                    {
                        biasx[ii] = 0.2 * biasx[ii] + 0.4 * (difx[ii] + biasx[ii]);
                        biasy[ii] = 0.2 * biasy[ii] + 0.4 * (dify[ii] + biasy[ii]);
                    }
                    numSuccess++;
                    NumFailed = 0;
                }
                else
                {
                    for (int ii = 0; ii < N; ii++)
                    {
                        Solzegonx[ii] = Solx[ii] - biasx[ii] - difx[ii];
                        Solzegony[ii] = Soly[ii] - biasy[ii] - dify[ii];
                    }
                    Node_treatment_localization(ref Solzegonx, ref Solzegony, N, rand);
                    fitness_Sol_zegon=compute_cost(Solzegonx, Solzegony, Anchorx, Anchory, range, range, rx, ry, noise_factor, gaussian);
                    if (fitness_Sol_zegon < fitness_Sol)
                    {
                        Solx = (double[])Solzegonx.Clone();
                        Soly = (double[])Solzegony.Clone();
                        fitness_Sol = fitness_Sol_zegon;
                        for (int ii = 0; ii < N; ii++)
                        {
                            biasx[ii] = biasx[ii] - 0.4 * (difx[ii] + biasx[ii]);
                            biasx[ii] = biasx[ii] - 0.4 * (difx[ii] + biasx[ii]);
                        }
                        numSuccess++;
                        NumFailed = 0;

                    }
                    else
                    {
                        NumFailed++;
                        numSuccess = 0;
                    }

                }

                if (numSuccess > 5)
                {
                    rho *= 2;
                    numSuccess = 0;
                }
                else if (NumFailed > 3)
                {
                    rho /= 2;
                    NumFailed = 0;
                }
            }

        }
        public LS[] Remove_Unchanged_particle(LS[] LS_S)
        {
            LS[] Result = new LS[LS_S.Length-1];
            for (int ii = 0; ii < Result.Length; ii++)
                Result[ii] = LS_S[ii + 1];
            return Result;
        }
        public void SolisWets(ref double[] Sol, ref double fitness_Sol, int N, ref  double[] bias, ref double rho, int MaxEval, Random rand, int OF,double[]xopt)
        {
        }
        public void SolisWets(ref double[] Sol, ref double fitness_Sol, int N, ref  double[] bias, ref double rho, int MaxEval, Random rand, int OF,double[]xopt,int dimension,int[]Perm,double[]r25,double[]r50,double[]r100,int[]S_Array,double[]W_Array,int overlap,int ub,int lb)
        {
            function_evaluation2013 class_new = new function_evaluation2013(ub,lb,dimension);
            int numSuccess = 0;
            int NumFailed = 0;
            double[] dif = new double[N];
            double[] Solprim = new double[N];
            double[] Solzegon = new double[N];
            double fitness_Sol_prim = 0;
            double fitness_Sol_zegon = 0;
            for (int numEval = 0; numEval < MaxEval; numEval++)
            {
                for (int ii = 0; ii < N; ii++)
                {
                    dif[ii] = Normal_Distribution(rand, 0, rho);
                    Solprim[ii] = Sol[ii] + bias[ii] + dif[ii];
                }
                Solprim = Node_treatment_new(Solprim, N, rand, 100, -100);
                fitness_Sol_prim = class_new.Choosing_Function(OF, Solprim, dimension, xopt, Perm, r25, r50, r100, S_Array, W_Array, overlap);
                if (fitness_Sol > fitness_Sol_prim)
                {
                    Sol = (double[])Solprim.Clone();
                    fitness_Sol = fitness_Sol_prim;
                    for (int ii = 0; ii < N; ii++)
                    {
                        bias[ii] = 0.2 * bias[ii] + 0.4 * (dif[ii] + bias[ii]);
                    }
                    numSuccess++;
                    NumFailed = 0;
                }
                else
                {
                    for (int ii = 0; ii < N; ii++)
                        Solzegon[ii] = Sol[ii] - bias[ii] - dif[ii];
                    Solzegon = Node_treatment(Solzegon, N, rand);
                    fitness_Sol_zegon = class_new.Choosing_Function(OF, Solprim, dimension, xopt, Perm, r25, r50, r100, S_Array, W_Array, overlap);
                    if (fitness_Sol_zegon < fitness_Sol)
                    {
                        Sol = (double[])Solzegon.Clone();
                        fitness_Sol = fitness_Sol_zegon;
                        for (int ii = 0; ii < N; ii++)
                        {
                            bias[ii] = bias[ii] - 0.4 * (dif[ii] + bias[ii]);
                        }
                        numSuccess++;
                        NumFailed = 0;

                    }
                    else
                    {
                        NumFailed++;
                        numSuccess = 0;
                    }

                }

                if (numSuccess > 5)
                {
                    rho *= 2;
                    numSuccess = 0;
                }
                else if (NumFailed > 3)
                {
                    rho /= 2;
                    NumFailed = 0;
                }
            }

        }
        public double[] Node_treatment(double[] Node, int N, Random rr)
        {
            double min = -1;
            double max = +1;
            for (int ii = 0; ii < N; ii++)
            {
                if (Node[ii] > max)
                {
                    Node[ii] = rr.NextDouble();
                    int temp = rr.Next(0, 2);
                    if (temp == 1)
                        Node[ii] = -Node[ii];
                }
                if (Node[ii] < min)
                {
                    Node[ii] = rr.NextDouble();
                    int temp = rr.Next(0, 2);
                    if (temp == 1)
                        Node[ii] = -Node[ii];
                }
            }
            return Node;
        }
        public double[] Node_treatment_new(double[] Node, int N, Random rr,double ma,double mi)
        {
            double min = mi;
            double max = ma;
            for (int ii = 0; ii < N; ii++)
            {
                if (Node[ii] > max)
                {
                    Node[ii] = rr.NextDouble();
                    int temp = rr.Next(0, 2);
                    if (temp == 1)
                        Node[ii] = -Node[ii];
                }
                if (Node[ii] < min)
                {
                    Node[ii] = rr.NextDouble();
                    int temp = rr.Next(0, 2);
                    if (temp == 1)
                        Node[ii] = -Node[ii];
                }
            }
            return Node;
        }
        public void Node_treatment(ref double[] x,ref double[] v, int N, Random rr)
        {
            double min = -1;
            double max = 1;
            for (int ii = 0; ii < N; ii++)
            {
                x[ii] = Math.Max(x[ii], min);
                x[ii] = Math.Min(x[ii], max);
                v[ii] = Math.Max(v[ii], min);
                v[ii] = Math.Min(v[ii], max);
            }
        }
        public void Node_treatment_localization(ref double[] x,ref double[] y, int N, Random rr)
        {
            double min = 0;
            double max = 1;
            for (int ii = 0; ii < N; ii++)
            {
                if (x[ii] > max)
                {
                    x[ii] = rr.NextDouble();
                }
                if (x[ii] < min)
                {
                    x[ii] = rr.NextDouble();
                }
                if (y[ii] > max)
                {
                    y[ii] = rr.NextDouble();
                }
                if (y[ii] < min)
                {
                    y[ii] = rr.NextDouble();
                }
            }
        }
        public double Min(double x1, double x2)
        {
            double min = 1000;
            if (x1 < min)
                min = x1;
            if (x2 < min)
                min = x2;


            return min;
        }
        public double Max(double x1, double x2)
        {
            double min = -1000;

                if (x1 > min)
                    min = x1;
                if (x2 > min)
                    min = x2;

            return min;
        }
        public int negative_assortative_mating_strategy(double[] x, LS[] Pop,Random rr)
        {
            double max=-1000;
            int index=0;
            int[] p = new int[3];
            p[0] = rr.Next(Pop.Length);
            p[1] = rr.Next(Pop.Length);
            p[2] = rr.Next(Pop.Length);
            for (int ii = 0; ii < p.Length; ii++)
            {
                double temp = 0;
                for (int jj = 0; jj < x.Length; jj++)
                    temp += Math.Abs(Pop[p[ii]].x[jj] - x[jj]);
                if(max<temp)
                {
                    max=temp;
                    index=p[ii];
                }
            }
            return index;
        }
        public double[] Hill_climbing_crossover(double[] x1, double[] x2, int nit, int Noff, int OF, int N,int[] Perm,double[,] Orthogonal_Matrix,int M,double[] o)
        {
            double[] p1 = (double[])x1.Clone();
            double[] p2 = (double[])x2.Clone();
            double P1_fit = 0; double p2_fit = 0;
            double[] Obest = new double[p1.Length];
            double[] result = new double[N];
            double Obest_fitness = 1000000;
            for (int i = 0; i < nit; i++)
            {
                bool infil = false;
                for (int jj = 0; jj < Noff; jj++)
                {
                    double[] temp = BLX_Crossover(p1, p2);
                    double temp_fit = compute_fitness(temp, OF, N, Perm, Orthogonal_Matrix, M, o);
                    if (Obest_fitness > temp_fit)
                    {
                        Obest_fitness = temp_fit;
                        Obest = (double[])temp.Clone();
                        infil = true;
                    }
                }
                if (infil)
                {
                    P1_fit = compute_fitness(p1, OF, N, Perm, Orthogonal_Matrix, M, o);
                    p2_fit = compute_fitness(p2, OF, N, Perm, Orthogonal_Matrix, M, o);
                    if (P1_fit < p2_fit)
                    {
                        //if (p2_fit > Obest_fitness)
                        {
                            p2 = (double[])Obest.Clone();
                        }
                    }
                    else
                    {
                        //if (P1_fit > Obest_fitness)
                        {
                            p1 = (double[])Obest.Clone();
                        }
                    }
                }
            }
            //P1_fit = simple_function(p1, OF, N);
            //p2_fit = simple_function(p2, OF, N);
            //if (P1_fit < p2_fit)
            //    result =(double[]) p1.Clone();
            //else
            //    result = (double[])p2.Clone();
            Random rr = new Random();
            if (rr.NextDouble() < 0.5)
                result = (double[])p1.Clone();
            else
                result = (double[])p2.Clone();
            return result;
        }
        public double diversity(LS[] Pop, int pop_no, int N)
        {
            double[] Ni = new double[pop_no];
            double diver = 0;
            for (int ii = 1; ii < pop_no - 1; ii++)
                for (int jj = ii + 1; jj < pop_no; jj++)
                    for (int kk = 0; kk < N; kk++)
                    {
                        diver += Math.Abs(Pop[ii].x[kk] - Pop[jj].x[kk]);
                    }

            //for (int kk = 0; kk < N; kk++)
            //{
            //    double min = 10000000; double max = -1000000;
            //    for (int ii = 1; ii < pop_no - 1; ii++)
            //    {
            //        if (max < Pop[ii].x[kk])
            //            max = Pop[ii].x[kk];
            //        if (min > Pop[ii].x[kk])
            //            min = Pop[ii].x[kk];
            //    }
            //    diver += Math.Abs(max-min);
            //}
            return diver/(N*pop_no*(pop_no-1));
        }
        public int positive_assortative_mating_strategy(double[] x, LS[] Pop, Random rr, int[] S)
        {
            double min = 1000;
            int index = 0;
            int[] p = new int[2];
            int temp1 = rr.Next(S.Length);
            int temp2 = rr.Next(S.Length);
            p[0] = S[temp1];
            p[1] = S[temp2];
            for (int ii = 0; ii < p.Length; ii++)
            {
                double temp = 0;
                for (int jj = 0; jj < x.Length; jj++)
                    temp += Math.Abs(Pop[p[ii]].x[jj] - x[jj]);
                if (min > temp)
                {
                    min = temp;
                    index = p[ii];
                }
            }
            return index;
        }
        public int negative_assortative_mating_strategy(double[] x, LS[] Pop, Random rr,int[]S)
        {
            double max = -1000;
            int index = 0;
            int[] p = new int[3];
            int temp1= rr.Next(S.Length);
            int temp2=rr.Next(S.Length);
            int temp3 = rr.Next(S.Length);
            p[0] = S[temp1];
            p[1] = S[temp2];
            p[2] = S[temp3];
            for (int ii = 0; ii < p.Length; ii++)
            {
                double temp = 0;
                for (int jj = 0; jj < x.Length; jj++)
                    temp += Math.Abs(Pop[p[ii]].x[jj] - x[jj]);
                if (max < temp)
                {
                    max = temp;
                    index = p[ii];
                }
            }
            return index;
        }
        public int negative_assortative_mating_strategy(double[] x,double[] y ,indi[] Pop, Random rr)
        {
            double max = -1000;
            int index = 0;
            int[] p = new int[3];
            p[0] = rr.Next(Pop.Length);
            p[1] = rr.Next(Pop.Length);
            p[2] = rr.Next(Pop.Length);
            for (int ii = 0; ii < p.Length; ii++)
            {
                double temp = 0;
                for (int jj = 0; jj < x.Length; jj++)
                    temp += Math.Pow(Math.Abs(Pop[p[ii]].x[jj] - x[jj]), 2) + Math.Pow(Math.Abs(Pop[p[ii]].y[jj] - y[jj]), 2);
                if (max < temp)
                {
                    max = temp;
                    index = p[ii];
                }
            }
            
            return index;
        }
        public int negative_assortative_mating_strategy2(double[] x, LS[] Pop, Random rr, int[] S)
        {
            double max = -1000;
            int index = 0;
            for (int ii = 0; ii < S.Length; ii++)
            {
                double temp = 0;
                for (int jj = 0; jj < x.Length; jj++)
                    temp += Math.Abs(Pop[S[ii]].x[jj] - x[jj]);
                if (max < temp)
                {
                    max = temp;
                    index = S[ii];
                }
            }
            return index;
        }
        public int negative_assortative_mating_strategy2(double[] x,double[]y ,indi[] Pop, Random rr, int[] S)
        {
            double max = -1000;
            int index = 0;
            for (int ii = 0; ii < S.Length; ii++)
            {
                double temp = 0;
                for (int jj = 0; jj < x.Length; jj++)
                    temp += Math.Pow(Math.Abs(Pop[S[ii]].x[jj] - x[jj]), 2) + Math.Pow(Math.Abs(Pop[S[ii]].y[jj] - y[jj]), 2);
                if (max < temp)
                {
                    max = temp;
                    index = S[ii];
                }
            }
            return index;
        }
        public double[] BLX_Crossover(double[] x1, double[] x2)
        {
            Random rr = new Random();
            double[] breed = new double[x2.Length];
            for (int ii = 0; ii < x1.Length; ii++)
            {
                double CMin = Min(x1[ii], x2[ii]);
                double CMax = Max(x1[ii], x2[ii]);
                double I = (CMax - CMin)*0.5;
                double A1 = CMin-I;
                if (A1 < CMin) 
	            A1 = CMin;
                double B1 = CMax+I;
                if (B1 > CMax) 
	                B1 = CMax;
                breed[ii] = A1 + rr.NextDouble() * (B1 - A1);
            }
                
            return breed;
        }
        public void phermone_decay(ref double[] x, double Gdc)
        {
        }
        public int Selection_Operators(int method,NF[]Pop,double min,int pop_no,Random Rand)
        {
            int result=-1;
            if (method == 0)
                result = Select_with_Probability(Pop, min, pop_no, Rand);
            else if (method == 1)
                result = ranking_selection(Pop, Rand, pop_no);
            else if (method == 2)
                result = Tournament(Pop, pop_no, Rand);
            //else if (method == 3)
            //    result = negative_assortative_mating(Pop, Rand);
            return result;
        }
        public int Tournament(NF[] Pop, int pop_no, Random Rand)
        {
            int r1 = Rand.Next(pop_no);
            int r2 = Rand.Next(pop_no);
            int index = 0;
            double max = -100;
            if (r1 > r2)
            {
                for(int ii=r2;ii<r1;ii++)
                    if (max < Pop[ii].F)
                    {
                        max = Pop[ii].F;
                        index = ii;
                    }
            }
            else
            {
                for (int ii = r1; ii < r2; ii++)
                    if (max < Pop[ii].F)
                    {
                        max = Pop[ii].F;
                        index = ii;
                    }
            }
            return index;
        }
        public int ranking_selection(NF[]Pop, Random rr,int pop_no)
        {
            double[] temp = new double[pop_no];
            int[] rank = new int[pop_no];
            int result = -1;
            for (int ii = 0; ii < pop_no; ii++)
            {
                rank[ii] = ii;
                temp[ii] = Pop[ii].F;
            }

            for (int ii = 0; ii < temp.Length - 1; ii++)
                for (int jj = ii + 1; jj < temp.Length; jj++)
                {
                    if (temp[ii] < temp[jj])
                    {
                        double temps = temp[ii];
                        temp[ii] = temp[jj];
                        temp[jj] = temps;
                        int rtemp = rank[ii];
                        rank[ii] = rank[jj];
                        rank[jj] = rtemp;
                    }
                }
            double p = 1;
            for (int ii = 0; ii < rank.Length; ii++)
            {
                temp[ii] = Math.Pow((1 - ((Convert.ToDouble(ii) + 1) / (rank.Length + 1))), p);
            }
            for (int ii = 0; ii < temp.Length; ii++)
                if (temp[ii] > rr.NextDouble())
                {
                    result = rank[ii];
                    break;
                }
            if (result == -1)
                result = rank[rank.Length - 1];
            return result;
        }
        public int negative_assortative_mating(NF[] Pop, Random rr)
        {
            double max = -1000;
            int index = 0;
            int rand_x = rr.Next(Pop.Length);
            double[] x =(double[])Pop[rand_x].x.Clone();
            int[] p = new int[3];
            p[0] = rr.Next(Pop.Length);
            p[1] = rr.Next(Pop.Length);
            p[2] = rr.Next(Pop.Length);
            for (int ii = 0; ii < p.Length; ii++)
            {
                double temp = 0;
                for (int jj = 0; jj < x.Length; jj++)
                    temp += Math.Abs(Pop[p[ii]].x[jj] - x[jj]);
                if (max < temp)
                {
                    max = temp;
                    index = ii;
                }
            }
            return index;
        }
        public NF Crossover_Operators(NF indi1, NF indi2,int N,int OF, Random Rand,int method)
        {
            NF result = new NF();
            if (method == 0)
                result = Cross_Over_one_point(indi1,indi2, N, OF, Rand);
            else if (method == 1)
                result = Cross_Over_two_point(indi1, indi2, N, OF, Rand);
            else if (method == 2)
                result = mask(indi1, indi2, N, OF, Rand);
            return result;
        }
        public double update_operator_SRE(NF[]Pop,double[]OldCost,double thresold)
        {
            double improvement_rate=0;
           for(int ii=0;ii<Pop.Length;ii++)
               improvement_rate+=(Pop[ii].F-OldCost[ii]);
            double mean=improvement_rate/Pop.Length;
            return mean;
        }
        public double[] Laplacian_recombination(LS par1, LS par2, int D, Random rand, double a, double b)
        {
            double[] result = new double[D];
            double[] y1 = new double[D];
            double[] y2 = new double[D];
            double u = rand.NextDouble(); double Betha = 0;
            if (u <= 0.5)
                Betha = a - b * Math.Log(u, Math.E);
            else
                Betha = a + b * Math.Log(u, Math.E);
            for (int ii = 0; ii < D; ii++)
            {
                y1[ii] = par1.x[ii] + Betha * Math.Abs(par1.x[ii] - par2.x[ii]);
                y2[ii] = par2.x[ii] + Betha * Math.Abs(par1.x[ii] - par2.x[ii]);
            }
            if (rand.NextDouble() < 0.5)
                result = (double[])y1.Clone();
            else
                result = (double[])y2.Clone();
            return result;

        }
        public void Subgrouping_SolisWets(ref double[] Sol, ref double fitness_Sol, int N, ref  double[] bias, ref double rho, int MaxEval, Random rand, int OF, int[] Perm, double[,] Orthogonal_Matrix, int M, double[] o)
        {
            int numSuccess = 0;
            int NumFailed = 0;
            int maxVars = 50;
            double[] dif = new double[N];
            double[] Solprim = new double[N];
            double[] Solzegon = new double[N];
            double fitness_Sol_prim = 0;
            double fitness_Sol_zegon = 0;
            bool[] change = new bool[N];
            int max_eval_subset = MaxEval / 10+1;
            for (int numEval = 0; numEval < MaxEval; numEval++)
            {
                if (Convert.ToInt32(numEval % max_eval_subset) == 0)
                {
                    change = new bool[N];
                    dif = new double[N];
                    set_Sub_set(ref change, N, maxVars);
                }
                for (int ii = 0; ii < N; ii++)
                {
                    if (change[ii] == true)
                    {
                        dif[ii] = Normal_Distribution(rand, 0, rho);
                    }
                    Solprim[ii] = Sol[ii] + bias[ii] + dif[ii];

                }
                Solprim = Node_treatment(Solprim, N, rand);
                if (OF < 16 || OF == 22)
                    fitness_Sol_prim = compute_fitness(Solprim, OF, N);
                else
                    fitness_Sol_prim = compute_fitness(Solprim, OF, N, Perm, Orthogonal_Matrix, M, o);
                if (fitness_Sol < fitness_Sol_prim)
                {
                    Sol = (double[])Solprim.Clone();
                    fitness_Sol = fitness_Sol_prim;
                    for (int ii = 0; ii < N; ii++)
                    {
                        bias[ii] = 0.2 * bias[ii] + 0.4 * (dif[ii] + bias[ii]);
                    }
                    numSuccess++;
                    NumFailed = 0;
                }
                else
                {
                    for (int ii = 0; ii < N; ii++)
                        Solzegon[ii] = Sol[ii] - bias[ii] - dif[ii];
                    Solzegon = Node_treatment(Solzegon, N, rand);
                    if (OF < 16 || OF == 22)
                        fitness_Sol_zegon = compute_fitness(Solzegon, OF, N);
                    else
                        fitness_Sol_zegon = compute_fitness(Solzegon, OF, N, Perm, Orthogonal_Matrix, M, o);
                    if (fitness_Sol_zegon > fitness_Sol)
                    {
                        Sol = (double[])Solzegon.Clone();
                        fitness_Sol = fitness_Sol_zegon;
                        for (int ii = 0; ii < N; ii++)
                        {
                            bias[ii] = bias[ii] - 0.4 * (dif[ii] + bias[ii]);
                        }
                        numSuccess++;
                        NumFailed = 0;

                    }
                    else
                    {
                        NumFailed++;
                        numSuccess = 0;
                    }

                }

                if (numSuccess > 5)
                {
                    rho *= 2;
                    numSuccess = 0;
                }
                else if (NumFailed > 3)
                {
                    rho /= 2;
                    NumFailed = 0;
                }
            }

        }
        public void set_Sub_set(ref bool[] change, int N, int maxVar)
        {
            Random r = new Random();
            int iniVar = r.Next(N - 1);
            int numVar = Convert.ToInt32(0.2 * N);
            if (numVar > maxVar)
                numVar = maxVar;
            for (int CurrentVar = iniVar; CurrentVar < (iniVar + numVar); CurrentVar++)
            {
                change[(CurrentVar % N)] = true;
            }
        }
        public void GetSpeciesCounts(ref BBO[] Pop)
        {
            for (int ii = 0; ii < Pop.Length; ii++)
            {
                Pop[ii].Species = Pop.Length - ii;
            }
        }
        public LS GetLambdaMu(BBO[]Pop,int I, int E)
        {
            LS lam_mu = new LS();
            double[] lambda = new double[Pop.Length];
            double[] mu = new double[Pop.Length];
            for (int ii = 0; ii < Pop.Length; ii++)
            {
                lambda[ii] = I * (1 - Pop[ii].Species / Pop.Length);
                mu[ii] = E * Pop[ii].Species / Pop.Length;
            }
            lam_mu.bias = (double[])lambda.Clone();
            lam_mu.x = (double[])mu.Clone();
            return lam_mu;
        }
        public bool similarity(double[] x1, double[] x2)
        {
            for (int ii = 0; ii < x1.Length; ii++)
                if (x1[ii] != x2[ii])
                    return false;
            return true;
        }
        public double mean_simple(ArrayList A)
        {
            double sum = 0; double sum2 = 0;
            for (int ii = 0; ii < A.Count; ii++)
            {
                sum += Convert.ToDouble(A[ii]);
                sum2 += Math.Pow(Convert.ToDouble(A[ii]), 2);
            }
            return sum2 / sum;
        }
        public double mean_Lehmer(ArrayList A)
        {
            double sum = 0;
            for (int ii = 0; ii < A.Count; ii++)
            {
                sum += Convert.ToDouble(A[ii]);
            }
            return sum / A.Count;
        }
        public void rosenbrock(int budget, double alpha, double beta, int problemDimension, ref double fBest, ref double[] best, Random Rand, int OF, int[] Perm, double[,] Orthogonal_Matrix, int M, double[] o, ref int it)
        {
            double eps = 1E-5; //10e-5

            double[] x = new double[problemDimension];
            double[,] xi = new double[problemDimension, problemDimension];
            for (int ii = 0; ii < problemDimension; ii++)
                xi[ii, ii] = 1;
            double[,] A = new double[problemDimension, problemDimension];
            double[] lambda = new double[problemDimension];
            double[] xCurrent = new double[problemDimension];
            double[] t = new double[problemDimension];
            double[] xk = new double[problemDimension];
            double[] d = new double[problemDimension];

            double yFirstFirst = fBest;
            double yFirst = yFirstFirst;
            double yCurrent;
            for (int n = 0; n < problemDimension; n++)
            {
                d[n] = 0.1;
                xk[n] = best[n];
                x[n] = best[n];
            }

            bool restart = true;
            bool firstExternalWhile = true;
            bool firstInternalWhile;
            double mini; double div;

            while (firstExternalWhile || (restart && it < budget))
            {
                fBest = yFirstFirst;
                firstInternalWhile = true;
                while (firstInternalWhile || (fBest > yFirst))
                {
                    yFirst = fBest;
                    for (int j = 0; j < problemDimension; j++)
                    {

                        for (int n = 0; n < problemDimension; n++)
                            xCurrent[n] = xk[n] + d[j] * xi[j, n];
                        xCurrent = Node_treatment(xCurrent, problemDimension, Rand);
                        if (OF < 16 || OF == 22)
                            yCurrent = compute_fitness(xCurrent, OF, problemDimension);
                        else
                            yCurrent = compute_fitness(xCurrent, OF, problemDimension, Perm, Orthogonal_Matrix, M, o);
                        it++;
                        if (it >= budget)
                            break;

                        if (yCurrent > fBest)
                        {
                            lambda[j] = lambda[j] + d[j];
                            d[j] = alpha * d[j];
                            fBest = yCurrent;
                            for (int n = 0; n < problemDimension; n++)
                            {
                                xk[n] = xCurrent[n];
                                best[n] = xCurrent[n];
                            }
                        }
                        else
                            d[j] = -beta * d[j];
                    }

                    if (it >= budget)
                        break;
                    firstInternalWhile = false;
                }

                mini = Min_all(d);
                restart = mini > eps;

                if (fBest > yFirstFirst)
                {
                    mini = Min_all(Subtract(xk, x));
                    restart = restart || (mini > eps);
                    if (restart)
                    {
                        for (int n = 0; n < problemDimension; n++)
                            A[problemDimension - 1, n] = lambda[n] * xi[problemDimension - 1, n];
                        t[problemDimension - 1] = lambda[problemDimension - 1] * lambda[problemDimension - 1];
                        for (int k = problemDimension - 2; k > 0; k--)
                        {
                            for (int n = 0; n < problemDimension; n++)
                                A[k, n] = A[k + 1, n] + lambda[n] * xi[k, n];
                            t[k] = t[k + 1] + lambda[k] * lambda[k];
                        }

                        for (int k = problemDimension - 1; k > 1; k--)
                        {
                            div = Math.Sqrt(t[k - 1] * t[k]);
                            if (div != 0)
                                for (int n = 0; n < problemDimension; n++)
                                    xi[k, n] = (lambda[k - 1] * A[k, n] - xi[k - 1, n] * t[k]) / div;
                        }
                        div = Math.Sqrt(t[0]);
                        for (int n = 0; n < problemDimension; n++)
                        {
                            if (div != 0)
                                xi[0, n] = A[0, n] / div;
                            x[n] = xk[n];
                            lambda[n] = 0;
                            d[n] = 0.1;
                        }
                        yFirstFirst = fBest;
                    }
                }

                firstExternalWhile = false;
            }
        }
        public void rosenbrock_localization(int budget, double alpha, double beta, int problemDimension, ref double fBest, ref double[] bestx, ref double[] besty, Random Rand, double range, double[] rx, double[] ry, double[] gaussian, double noise_factor, double[] Anchorx, double[] Anchory, ref int it)
        {
            double eps = 1E-5; //10e-5

            double[] x = new double[problemDimension];
            double[] y = new double[problemDimension];
            double[,] xi = new double[problemDimension, problemDimension];
            double[,] yi = new double[problemDimension, problemDimension];
            for (int ii = 0; ii < problemDimension; ii++)
            {
                xi[ii, ii] = 1;
            }
            double[,] A = new double[problemDimension, problemDimension];
            double[] lambda = new double[problemDimension];
            double[] xCurrent = new double[problemDimension];
            double[] yCurrent = new double[problemDimension];
            double[] t = new double[problemDimension];
            double[] xk = new double[problemDimension];
            double[] yk = new double[problemDimension];
            double[] d = new double[problemDimension];

            double yFirstFirst = fBest;
            double yFirst = yFirstFirst;
            double fCurrent;
            for (int n = 0; n < problemDimension; n++)
            {
                d[n] = 0.1;
                xk[n] = bestx[n];
                yk[n] = besty[n];
                x[n] = bestx[n];
                y[n] = besty[n];
            }

            bool restart = true;
            bool firstExternalWhile = true;
            bool firstInternalWhile;
            double mini; double div;

            while (firstExternalWhile || (restart && it < budget))
            {
                fBest = yFirstFirst;
                firstInternalWhile = true;
                while (firstInternalWhile || (fBest > yFirst))
                {
                    yFirst = fBest;
                    for (int j = 0; j < problemDimension; j++)
                    {

                        for (int n = 0; n < problemDimension; n++)
                        {
                            xCurrent[n] = xk[n] + d[j] * xi[j, n];
                            yCurrent[n] = yk[n] + d[j] * yi[j, n];
                        }
                        Node_treatment_localization(ref xCurrent,ref yCurrent, problemDimension, Rand);
                        fCurrent = compute_cost(xCurrent, yCurrent, Anchorx, Anchory, range, range, rx, ry, noise_factor, gaussian);
                        it++;
                        if (it >= budget)
                            break;

                        if (fCurrent > fBest)
                        {
                            lambda[j] = lambda[j] + d[j];
                            d[j] = alpha * d[j];
                            fBest = fCurrent;
                            for (int n = 0; n < problemDimension; n++)
                            {
                                xk[n] = xCurrent[n];
                                yk[n] = yCurrent[n];
                                bestx[n] = xCurrent[n];
                                besty[n] = yCurrent[n];
                            }
                        }
                        else
                            d[j] = -beta * d[j];
                    }

                    if (it >= budget)
                        break;
                    firstInternalWhile = false;
                }

                mini = Min_all(d);
                restart = mini > eps;

                if (fBest > yFirstFirst)
                {
                    mini = Min_all(Subtract(xk, x));
                    restart = restart || (mini > eps);
                    if (restart)
                    {
                        for (int n = 0; n < problemDimension; n++)
                        {
                            A[problemDimension - 1, n] = lambda[n] * xi[problemDimension - 1, n];
                        }
                        t[problemDimension - 1] = lambda[problemDimension - 1] * lambda[problemDimension - 1];
                        for (int k = problemDimension - 2; k > 0; k--)
                        {
                            for (int n = 0; n < problemDimension; n++)
                                A[k, n] = A[k + 1, n] + lambda[n] * xi[k, n];
                            t[k] = t[k + 1] + lambda[k] * lambda[k];
                        }

                        for (int k = problemDimension - 1; k > 1; k--)
                        {
                            div = Math.Sqrt(t[k - 1] * t[k]);
                            if (div != 0)
                                for (int n = 0; n < problemDimension; n++)
                                    xi[k, n] = (lambda[k - 1] * A[k, n] - xi[k - 1, n] * t[k]) / div;
                        }
                        div = Math.Sqrt(t[0]);
                        for (int n = 0; n < problemDimension; n++)
                        {
                            if (div != 0)
                                xi[0, n] = A[0, n] / div;
                            x[n] = xk[n];
                            lambda[n] = 0;
                            d[n] = 0.1;
                        }
                        yFirstFirst = fBest;
                    }
                }

                firstExternalWhile = false;
            }
        }
        public double[] Subtract(double[] x, double[] xk)
        {
            double[] result = new double[x.Length];
            for (int ii = 0; ii < x.Length; ii++)
                result[ii] = x[ii] - xk[ii];
            return result;
        }
        public double Min_all(double[] h)
        {
            double min = 1000;
            double[] h_temp = (double[])h.Clone();
            for (int ii = 0; ii < h.Length; ii++)
            {
                h_temp[ii] = Math.Abs(h_temp[ii]);
                if (min > h_temp[ii])
                    min = h_temp[ii];
            }
            return min;
        }
        public double[] Construction_space(int Swarm_size,int number_of_swarm,Swarms[] SW, int index_particle,int N)
        {
            double[] x=new double[N];int k=0;
            for(int ii=0; ii<number_of_swarm;ii++)
                for (int jj = 0; jj < Swarm_size; jj++)
                {
                    x[k++] = SW[ii].particles[index_particle].x[jj];
                }
            return x;
        }
        public double[] Concatenation(double[] best,int N ,int position, int swarm_size, double[] x)
        {
            double[] result =(double[])best.Clone();
            int index = 0;
            for (int ii = position; ii <position+ swarm_size; ii++)
            {
                result[ii] = x[index++];
            }
            return result;
        }
        public double Cauchy_mu(Random rr)
        {
            return Math.Tan(Math.PI * (rr.NextDouble() - 0.5));
        }
        public indi construct_Non_Anchor_node(int NNA, Random rr, int Area_Size)
        {
            indi A = new indi();
            double[] Ax = new double[NNA];
            double[] Ay = new double[NNA];
            for (int ii = 0; ii < NNA; ii++)
            {
                Ax[ii] = rr.NextDouble() * Area_Size;
                Ay[ii] = rr.NextDouble() * Area_Size;
            }
            A.x = (double[])Ax.Clone();
            A.y = (double[])Ay.Clone();
            return A;
        }
        public indi Initialization(int No_nonAnch, Random RGN)
        {
            double[] x = new double[No_nonAnch];
            double[] y = new double[No_nonAnch];
            indi all = new indi();
            for (int ii = 0; ii < No_nonAnch; ii++)
            {
                x[ii] = RGN.NextDouble();
                y[ii] = RGN.NextDouble();
            }
            all.x = (double[])x.Clone();
            all.y = (double[])y.Clone();
            return all;
        }
        public double compute_cost(double[] Esx, double[] Esy, double[] Ax, double[] Ay, double rNA, double rA, double[] Rex, double[] Rey, double Noise_Factor, double[] gaussain)
        {
            double sum = 0;
            for (int jj = 0; jj < Esx.Length; jj++)
            {
                for (int ii = 0; ii < Ax.Length; ii++)
                {
                    double temp = compute_distance(Ax[ii], Ay[ii], Rex[jj], Rey[jj]);
                    if (temp <= rA)
                    {
                        sum += Math.Pow(compute_distance(Ax[ii], Ay[ii], Esx[jj], Esy[jj]) - measured_distance(compute_distance(Ax[ii], Ay[ii], Rex[jj], Rey[jj]), gaussain[jj], Noise_Factor), 2);
                    }
                }

                for (int ii = 0; ii < Esx.Length; ii++)
                    if (Esx[jj] != Esx[ii])
                    {
                        double temp = compute_distance(Rex[ii], Rey[ii], Rex[jj], Rey[jj]);
                        if (temp <= rNA)
                        {
                            sum += Math.Pow(compute_distance(Esx[ii], Esy[ii], Esx[jj], Esy[jj]) - measured_distance(compute_distance(Rex[ii], Rey[ii], Rex[jj], Rey[jj]), gaussain[jj], Noise_Factor), 2);
                        }
                    }
            }
            return sum;
        }
        public double compute_distance(double x1, double y1, double x2, double y2)
        {
            return Math.Sqrt(Math.Pow(x1 - x2, 2) + Math.Pow(y1 - y2, 2));
        }
        public double measured_distance(double distance, double guassian_variable, double noise_factor)
        {
            return distance * (1 + guassian_variable * noise_factor);
        }
        public double Calc_error(double[] Estx, double[] Esty, double[] Realx, double[] Realy, double range)
        {
            double sum = 0;
            for (int ii = 0; ii < Estx.Length; ii++)
                sum += Math.Pow(Estx[ii] - Realx[ii], 2) + Math.Pow(Esty[ii] - Realy[ii], 2);
            return 100 * sum / (Math.Pow(range, 2) * Realx.Length);

        }
        public void Subgrouping_SolisWets_localization(ref double[] Solx, ref double[] Soly, ref double fitness_Sol, int N, ref  double[] biasx, ref double[] biasy, ref double rho, int MaxEval, Random rand, double range, double[] rx, double[] ry, double[] gaussian, double noise_factor, double[] Anchorx, double[] Anchory)
        {
            int numSuccess = 0;
            int NumFailed = 0;
            int maxVars = 50;
            double[] difx = new double[N];
            double[] dify = new double[N];
            double[] Solprimx = new double[N];
            double[] Solprimy = new double[N];
            double[] Solzegonx = new double[N];
            double[] Solzegony = new double[N];
            double fitness_Sol_prim = 0;
            double fitness_Sol_zegon = 0;
            bool[] change = new bool[N];
            int max_eval_subset = MaxEval / 10 + 1;
            for (int numEval = 0; numEval < MaxEval; numEval++)
            {
                if (Convert.ToInt32(numEval % max_eval_subset) == 0)
                {
                    change = new bool[N];
                    difx = new double[N];
                    dify = new double[N];
                    set_Sub_set(ref change, N, maxVars);
                }
                for (int ii = 0; ii < N; ii++)
                {
                    if (change[ii] == true)
                    {
                        difx[ii] = Normal_Distribution(rand, 0, rho);
                        dify[ii] = Normal_Distribution(rand, 0, rho);
                    }
                    Solprimx[ii] = Solx[ii] + biasx[ii] + difx[ii];
                    Solprimy[ii] = Soly[ii] + biasy[ii] + dify[ii];

                }
                Node_treatment_localization(ref Solprimx, ref Solprimy, N, rand);
                fitness_Sol_prim = compute_cost(Solprimx, Solprimy, Anchorx, Anchory, range, range, rx, ry, noise_factor, gaussian);
                if (fitness_Sol > fitness_Sol_prim)
                {
                    Solx = (double[])Solprimx.Clone();
                    Soly = (double[])Solprimy.Clone();
                    fitness_Sol = fitness_Sol_prim;
                    for (int ii = 0; ii < N; ii++)
                    {
                        biasx[ii] = 0.2 * biasx[ii] + 0.4 * (difx[ii] + biasx[ii]);
                        biasy[ii] = 0.2 * biasy[ii] + 0.4 * (dify[ii] + biasy[ii]);
                    }
                    numSuccess++;
                    NumFailed = 0;
                }
                else
                {
                    for (int ii = 0; ii < N; ii++)
                    {
                        Solzegonx[ii] = Solx[ii] - biasx[ii] - difx[ii];
                        Solzegony[ii] = Soly[ii] - biasy[ii] - dify[ii];
                    }
                    Node_treatment_localization(ref Solzegonx, ref Solzegony, N, rand);
                    fitness_Sol_zegon = compute_cost(Solzegonx, Solzegony, Anchorx, Anchory, range, range, rx, ry, noise_factor, gaussian);
                    if (fitness_Sol_zegon > fitness_Sol)
                    {
                        Solx = (double[])Solzegonx.Clone();
                        Soly = (double[])Solzegony.Clone();
                        fitness_Sol = fitness_Sol_zegon;
                        for (int ii = 0; ii < N; ii++)
                        {
                            biasx[ii] = biasx[ii] - 0.4 * (difx[ii] + biasx[ii]);
                            biasx[ii] = biasx[ii] - 0.4 * (difx[ii] + biasx[ii]);
                        }
                        numSuccess++;
                        NumFailed = 0;

                    }
                    else
                    {
                        NumFailed++;
                        numSuccess = 0;
                    }

                }

                if (numSuccess > 5)
                {
                    rho *= 2;
                    numSuccess = 0;
                }
                else if (NumFailed > 3)
                {
                    rho /= 2;
                    NumFailed = 0;
                }
            }
            
        }
        public void Uniformly_distribution(ref double[] Realx, ref double[] Realy, Random rr, int NNA, int NA, double problem_space)
        {
            double region_sizex = problem_space / 5;
            double region_sizey = problem_space / 4;
            int jj = 0;
            int kk = 0;
            int numberx=4;
            int numbery=5;
            int number = NNA / NA;
            for (int ii=0;kk<NNA ; )
            {
                Realx[kk] = ii * region_sizex + region_sizex * rr.NextDouble();
                Realy[kk++] = jj * region_sizey + region_sizey * rr.NextDouble();
                if ((kk % number)==0&&kk!=0)
                {
                    jj++;
                    if (jj == numberx)
                    {
                        jj = 0;
                        ii++;
                        if (ii == numbery)
                            break;
                    }
                }
                
            }


        }
        public void Uniformly_distribution_anchor_nodes(ref double[] Anchorx,ref double[] Anchory, Random rr, int NNA, int NA, double problem_space)
        {
            double region_sizex = problem_space / 5;
            double region_sizey = problem_space / 4;
            int jj = 0;
            int kk = 0;
            int numberx = 4;
            int numbery = 5;
            int number = NA / NA;
            for (int ii = 0;; )
            {
                Anchorx[kk] = ii * region_sizex + region_sizex * rr.NextDouble();
                Anchory[kk++] = jj * region_sizey + region_sizey * rr.NextDouble();
                if ((kk % number) == 0 && kk != 0)
                {
                    jj++;
                    if (jj == numberx)
                    {
                        jj = 0;
                        ii++;
                        if (ii == numbery)
                            break;
                    }
                }

            }

        }
        public double[] generating_jumping(indi[] Nodes, int pop, int N, int index)
        {
            double min = 1000;
            double max = -100;
            double[] ONodes = new double[N];
            for (int jj = 0; jj < N; jj++)
            {
                for (int ii = 0; ii < pop; ii++)
                {
                    if (min > Nodes[ii].x[jj])
                        min = Nodes[ii].x[jj];
                    if (max < Nodes[ii].x[jj])
                        max = Nodes[ii].x[jj];
                }
            }
            for (int jj = 0; jj < N; jj++)
                ONodes[jj] = min + max - Nodes[index].x[jj];
            return ONodes;
        }
        public double[] generating_jumping(NF[] Nodes, int pop, int N, int index)
        {
            double min = 1000;
            double max = -100;
            double[] ONodes = new double[N];
            for (int jj = 0; jj < N; jj++)
            {
                for (int ii = 0; ii < pop; ii++)
                {
                    if (min > Nodes[ii].x[jj])
                        min = Nodes[ii].x[jj];
                    if (max < Nodes[ii].x[jj])
                        max = Nodes[ii].x[jj];
                }
            }
            for (int jj = 0; jj < N; jj++)
                ONodes[jj] = min + max - Nodes[index].x[jj];
            return ONodes;
        }
    }
}
