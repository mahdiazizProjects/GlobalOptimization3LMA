using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Text;
using System.Windows.Forms;
using System.IO;
using System.Collections;
using System.Linq;
namespace OMOA_Namespace
{


    struct NF
    {
        public double[] x;
        public double F;
    }
    struct LS
    {
        public double[] x;
        public double F;
        public double[] bias;
        public double rho;
        public int visited;
        public bool local_search;
    }
    struct R
    {
        public double Cost;
        public int conf;
        public double value;
    }
    public partial class Form1 : Form
    {

        StreamWriter sw;
        public Form1()
        {
            InitializeComponent();
        }
        private void Solution_Treatment(int method,ref double[]x,Random rr)
        {
            if(method==1)
            {
                for(int ii=0;ii<x.Length;ii++)
                {
                    if(x[ii]<-1||x[ii]>1)
                    {
                        x[ii] = -1 +rr.NextDouble();
                    }
                }
            }
            else if(method==2)
            {
                for (int ii = 0; ii < x.Length; ii++)
                {
                    if (x[ii] < -1)
                    {
                        x[ii] = -1 +0.001* rr.NextDouble();
                    }
                    if (x[ii] >1)
                    {
                        x[ii] = 1 - 0.001 * rr.NextDouble();
                    }
                }
            }
            else if(method==3)
            {
                for (int ii = 0; ii < x.Length; ii++)
                {
                    if (x[ii] < -1)
                    {
                        x[ii] = 2 * x[ii] + 1;
                    }
                    if (x[ii] > 1)
                    {
                        x[ii] = 2 * x[ii] - 1;
                    }
                }
            }
        }
        private void Form1_Load(object sender, EventArgs e)
        {

            orgMOA.Checked = true;
            SA.Enabled = false;
            GA.Enabled = false;
            PSO.Enabled = false;
            MOA.Enabled = false;
            ODE.Enabled = false;
            Cm_ODE.Text = "0.2";
            Cmax_ODE.Text = "1";
            C_step_ODE.Text = "0.2";
            Fmin_ODE.Text = "0.2";
            Fmax_ODE.Text = "1";
            F_s_ODE.Text = "0.2";
            Jr_min_ODE.Text = "0.3";
            Jr_max_ODE.Text = "0.3";
            Jr_step.Text = "0.3";
            Mu1.Text = "0";
            Mu2.Text = "4";
            Mu_step.Text = "1";
            cr1.Text = "0.2";
            cr2.Text = "1";
            cr_step.Text = ".2";
            TExp_text.Text = ".005";
            Cexp_text.Text = "0.005";
            Exp.Text = "NO";
            CSRR_text.Text = "0.5";
            ISSR_combo.Text = "0";
            SHR.Text = "0.01";
            CE.Text = "1";
            D.Text = "1";
            cmin.Text = "0.9";
            cmax.Text = "0.9";
            c_step.Text = "0.1";
            cmin21.Text = "0.9";
            cmin22.Text = "0.9";
            c2_step.Text = "0.1";
            Wmin.Text = "0.9";
            Wmax.Text = "0.9";
            W_step.Text = "0.9";
            Ac.Text = "0";
            D.Text = "1";
            Inten.Text = ".6";
            Al.Text = ".2";
            Inten2.Text = "0.6";
            Al2.Text = "0.2";
            Inten_step.Text = ".2";
            Al_step.Text = ".2";
            NT.Text = "50";
            NI.Text = "1000";
            NP.Text = "50";
            NE.Text = "100";
            NE2.Text = "100";
            NE_step.Text = "1";
            comboBox2.Text = "Cellular";
            N1.Text = "20";
            N2.Text = "100";
            N_s.Text = "20";
            colrate1.Text = ".91";
            colrate2.Text = ".99";
            col_step.Text = "0.02";
            T1.Text = "100";
            T2.Text = "1000";
            T_step.Text = "900";
            LI.Text = "100";
            of1.Text = "1";
            of2.Text = "21";
            of_s.Text = "1";
            Jumping_rate.Text = "0.3";
            jumping_rate2.Text = "1";
            j_step.Text = "0.05";
        }

        private void button1_Click(object sender, EventArgs e)
        {
            #region Initialization
            //double ub = 10;
            //double lb = -10;
            //int Dimension = 1000;
            int boundary_method = 1;
            int Ntest = Convert.ToInt32(NT.Text);
            int NPop = Convert.ToInt32(NP.Text);
            int Iteration = Convert.ToInt32(NI.Text);
            int ND1 = Convert.ToInt32(NE.Text);
            int ND2 = Convert.ToInt32(NE2.Text);
            int dimension=-1;
            int lb=-1,ub=-1;
            int ND_step = Convert.ToInt32(NE_step.Text);
            int OF1 = Convert.ToInt32(of1.Text);
            int OF2 = Convert.ToInt32(of2.Text);
            int OF_s = Convert.ToInt32(of_s.Text);
            StreamReader SR2 = new StreamReader("t-table.txt");
            string file = SR2.ReadToEnd();
            string[] FF = file.Split('\n');
            double[,] t_table = new double[102, 7];
            int kkkk = 0;
            for (int ii = 0; ii < FF.Length; ii++)
            {
                string[] temp = FF[ii].Split(' ');
                for (int jj = 0; jj < temp.Length; jj++)
                    if (temp[jj] != "")
                        t_table[ii, kkkk++] = Convert.ToDouble(temp[jj]);
                kkkk = 0;
            }
            string structure = comboBox2.Text;
            string method = comboBox1.Text;
            Numerical_functions NFC = new Numerical_functions();
            #endregion
            if (Criterion.Text != "RLD")
            {
                #region Benchmark Functions
            for (int N = ND1; N <= ND2; N++)
                for (int OF = OF1; OF <= OF2; OF += OF_s)
                {
                    if (OF != 14)
                    {

                        Random rrr = new Random(N);
                        int M = 20;
                        int[] Perm = NFC.sq_rand_gen2(N, N, rrr);
                        double[,] Orthogonal_Matrix;
                        double[] o = new double[N];
                        Random RO = new Random(100);
                        o = NFC.node_constructing(N, RO);
                        if (OF > 22)
                        {
                            M = 50;
                            Orthogonal_Matrix = new double[M, M];
                            if (OF != 24 && OF < 27)
                            {
                                StreamReader SR = new StreamReader("m" + OF.ToString() + ".txt");
                                string all = SR.ReadToEnd();
                                string[] all_with_space = all.Split('\n', ',');

                                for (int ii = 0; ii < M; ii++)
                                    for (int jj = 0; jj < M; jj++)
                                    {
                                        Orthogonal_Matrix[ii, jj] = Convert.ToDouble(all_with_space[ii + jj]);
                                    }
                            }
                        }
                        else
                        {
                            if (OF != 16)
                            {
                                StreamReader SR = new StreamReader("orthog" + M.ToString() + ".txt");
                                string all = SR.ReadToEnd();
                                string[] all_with_space = all.Split('\n', ',');
                                Orthogonal_Matrix = new double[M, M];
                                for (int ii = 0; ii < M; ii++)
                                    for (int jj = 0; jj < M; jj++)
                                    {
                                        Orthogonal_Matrix[ii, jj] = Convert.ToDouble(all_with_space[ii + jj]);
                                    }
                            }
                            else
                            {
                                StreamReader SR = new StreamReader("orthog" + N.ToString() + ".txt");
                                string all = SR.ReadToEnd();
                                string[] all_with_space = all.Split('\n', ',');
                                Orthogonal_Matrix = new double[N, N];
                                for (int ii = 0; ii < N; ii++)
                                    for (int jj = 0; jj < N; jj++)
                                    {
                                        Orthogonal_Matrix[ii, jj] = Convert.ToDouble(all_with_space[ii + jj]);
                                    }
                            }
                        }
                        //string ss2 = "";
                        double[] costs = new double[Ntest];
                        NF[] Nodes_test = new NF[Ntest];
                        Random Rand = new Random();
                        double Fmax = -1; double Fmin = 1000;
                        int O=50000;
                        int[] Pops = new int[] { 100,500,1000 };
                        for (int ip = 0; ip < Pops.Length; ip++)
                        {
                            NPop = Pops[ip];
                            Iteration = O / NPop;
                            string result = "fit_N_" + N.ToString() + "_method_" + method + "_pop_" + NPop.ToString() + "_OF_" + OF.ToString() + ".txt";
                            #region GA
                            if (method == "Genetic")
                            {
                                //double[] mus = { 0.002, 0.003, 0.005, 0.01,0.5 };
                                //double[] Crs = { 0.2, 0.4, 0.6, 0.8, 1 };
                                //double[] mus = { 0.005, 0.5, 0.5, 0.5, 0.5, 0.5 };
                                //double[] Crs = { 1, 0.4, 1, 0.4, 0.2, 0.4 };
                                //for (int muindex = 0; muindex < mus.Length; muindex++)
                                //for (int cr_index = 0; cr_index < Crs.Length; cr_index++)
                                //{
                                double[] mus = { 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.005, 0.003, 0.003, 0.003, 0.005, 0.003, 0.003, 0, 0.003, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.005 };
                                double[] Crs = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.8, 0.4, 0.4, 1.0, 0.8, 0.2, 1.0 };
                                //double cr_rate = Crs[OF - 22];
                                //double mu_rate = mus[OF - 22];
                                double cr_rate = Crs[OF];
                                double mu_rate = mus[OF];
                                //double mu_rate = muindex[Convert.ToInt32(MM)];


                                //ss2 = "Results\\fit_N_" + N.ToString() + "_bound_" + boundary_method.ToString()  +"pop_" + NPop.ToString() + method + "_cr_" + cr_rate.ToString() + "_mu_rate_" + mu_rate.ToString() + "_test_" + Ntest.ToString() + "OF_" + OF.ToString() + ".txt";
                                #region file_existency
                                if (!File.Exists(result))
                                {
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        NF[] GA = new NF[NPop];
                                        NF Best_indi = new NF();
                                        Best_indi.F = Double.NegativeInfinity;
                                        Fmax = Double.NegativeInfinity; Fmin = 10000000000;
                                        for (int ii = 0; ii < NPop; ii++)
                                        {

                                            GA[ii].x = NFC.node_constructing(N, Rand);
                                            if (OF < 16 || OF == 22)
                                                GA[ii].F = NFC.compute_fitness(GA[ii].x, OF, N);
                                            else
                                                GA[ii].F = NFC.compute_fitness(GA[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                            if (Fmax < GA[ii].F)
                                                Fmax = GA[ii].F;
                                        }
                                        for (int it = 0; it < Iteration; it++)
                                        {
                                            Fmax = Double.NegativeInfinity; Fmin = 100000000;
                                            int max_index = -1;
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                if (OF < 16 || OF == 22)
                                                    GA[ii].F = NFC.compute_fitness(GA[ii].x, OF, N);
                                                else
                                                    GA[ii].F = NFC.compute_fitness(GA[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                if (Fmax < GA[ii].F)
                                                {
                                                    Fmax = GA[ii].F;
                                                    max_index = ii;
                                                }
                                                if (Fmin > GA[ii].F)
                                                    Fmin = GA[ii].F;
                                            }
                                            if (Fmin == Fmax)
                                            {
                                                for (int ii = 0; ii < NPop; ii++)
                                                {
                                                    GA[ii].x = NFC.node_constructing(N, Rand);
                                                    if (OF < 16 || OF == 22)
                                                        GA[ii].F = NFC.compute_fitness(GA[ii].x, OF, N);
                                                    else
                                                        GA[ii].F = NFC.compute_fitness(GA[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                    if (Fmax < GA[ii].F)
                                                        Fmax = GA[ii].F;
                                                    if (Fmin > GA[ii].F)
                                                        Fmin = GA[ii].F;
                                                }
                                            }
                                            if (Fmax > Best_indi.F)
                                            {
                                                Best_indi.F = Fmax;
                                                Best_indi.x = (double[])GA[max_index].x.Clone();
                                            }

                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                if (cr_rate > Rand.NextDouble())
                                                {
                                                    int index = NFC.Select_with_Probability(GA, Fmin, NPop, Rand);
                                                    int index2 = NFC.Select_with_Probability(GA, Fmin, NPop, Rand);
                                                    NF childs = NFC.Cross_Over_one_point(GA[index], GA[index2], N, OF, Rand);
                                                    GA[ii].x = (double[])childs.x.Clone();
                                                }
                                                NFC.mutation(ref GA[ii], mu_rate, N, Rand);
                                                Solution_Treatment(boundary_method, ref GA[ii].x, Rand);
                                            }
                                        }
                                        costs[tt] = Best_indi.F;
                                        Nodes_test[tt].x = (double[])Best_indi.x.Clone();
                                        Range.Text = mu_rate.ToString() + " " + cr_rate.ToString() + " " + OF.ToString() + " " + tt.ToString();
                                        Range.Refresh();
                                    }


                                    sw = new StreamWriter(result);
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        sw.Write(costs[tt].ToString() + " ");
                                    }
                                    sw.Close();
                                    //}
                                }
                                #endregion



                            }
                            #endregion
                            #region MOA
                            if (method == "MOA")
                            {
                                double alpha = 0;
                                double intensity = 0;
                                double short_range = 0.01;
                                int acceleration = Convert.ToInt32(Ac.Text);
                                int Distance = Convert.ToInt32(D.Text);
                                double[] Magnet = new double[NPop];
                                double[] Mass = new double[NPop];
                                double al1 = Convert.ToDouble(Al.Text);
                                double al2 = Convert.ToDouble(Al2.Text);
                                double al_step = Convert.ToDouble(Al_step.Text);
                                double in1 = Convert.ToDouble(Inten.Text);
                                double in2 = Convert.ToDouble(Inten2.Text);
                                double in_step = Convert.ToDouble(Inten_step.Text);
                                //double[] Intinsities = { 0.00001 };
                                //double[] mediator = { 0.00001, 0.001, 0.01, 1, 10 };
                                //bool change = false;
                                //double[] Intinsities = { 0.001, 0.01, 0.2, 0.6, 1 };
                                //double[] Alphas = { 0.001, 0.01, 0.2, 0.6, 1 };
                                //double[] Intinsities = { 0.01, 0.2, 0.001, 0.01, 0.2, 0.001 };
                                //double[] Alphas = { 1, 0.001, 0.2, 0.2, 0.001, 0.2 };
                                double[] Intinsities = { 0.0010, 0.0100, 0.0100, 0.0100, 0.0100, 0.0100, 10.0000, 0.0100, 0.0010, 0.0010, 0.0100, 0.0010, 0.0100, 0, 0.0100, 0.0010, 0.0010, 0.0100, 0.0010, 0.0010, 0.0100, 0.0010 };
                                double[] Alphas = { 0.0010, 1.0000, 0.0100, 1.0000, 1.0000, 1.0000, 0.0010, 1.0000, 1.0000, 1.0000, 0.0100, 1.0000, 1.0000, 0, 1.0000, 1.0000, 1.0000, 1.0000, 0.0010, 1.0000, 1.0000, 1.0000 };
                                //for (int alpha_index = 0; alpha_index < 5; alpha_index++)
                                //for (int iiii = 0; iiii < 5; iiii++)
                                //{
                                //if (iiii == 1)
                                //    change = true;

                                //for (int intensity_index = 0; intensity_index < 5; intensity_index++)
                                //{
                                //    if (change == false)
                                //    {
                                //alpha = Intinsities[OF - 22];
                                //intensity = Alphas[OF - 22];
                                alpha = Intinsities[OF];
                                intensity = Alphas[OF];
                                //}
                                //else
                                //{
                                //    intensity = Intinsities[0];
                                //    alpha = mediator[intensity_index];
                                //}
                                //double intensity = Intinsities[OF - 1];
                                //double alpha = Alphas[OF - 1];

                                //ss2 = "Results\\fit_N_" + N.ToString() + "_bound_" + boundary_method.ToString()  +"_Pop_" + NPop.ToString() + "_intensity" + intensity.ToString() + "alpha_" + alpha.ToString() + "_Dis_" + method + Distance.ToString() + "_Acc_" + Ac.Text + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_struct_" + structure + "_OF" + OF.ToString() + ".txt";
                                if (!File.Exists(result))
                                {
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        NF[] Nodes = new NF[NPop];
                                        NF[] Forces = new NF[NPop];
                                        NF[] a = new NF[NPop];
                                        NF[] v = new NF[NPop];
                                        NF Best_indi = new NF();
                                        double[] distances = new double[NPop];
                                        Best_indi.F = -1E280;
                                        for (int ii = 0; ii < NPop; ii++)
                                        {
                                            Forces[ii].x = new double[N];
                                            a[ii].x = new double[N];
                                            v[ii].x = new double[N];
                                        }
                                        for (int ii = 0; ii < NPop; ii++)
                                        {
                                            Nodes[ii].x = NFC.node_constructing(N, Rand);
                                        }
                                        for (int it = 0; it < Iteration; it++)
                                        {
                                            Fmax = -1E280; Fmin = 100000000;
                                            int min_index = 0; int max_index = 0;
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                if (OF < 16 || OF == 22)
                                                    Nodes[ii].F = NFC.compute_fitness(Nodes[ii].x, OF, N);
                                                else
                                                    Nodes[ii].F = NFC.compute_fitness(Nodes[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                if (Fmax < Nodes[ii].F)
                                                {
                                                    Fmax = Nodes[ii].F;
                                                    max_index = ii;
                                                }
                                                if (Fmin > Nodes[ii].F)
                                                {
                                                    Fmin = Nodes[ii].F;
                                                    min_index = ii;
                                                }
                                            }

                                            if (Best_indi.F < Fmax)
                                            {
                                                Best_indi.F = Fmax;
                                                Best_indi.x = (double[])Nodes[max_index].x.Clone();
                                            }
                                            double range = Fmax - Fmin;
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                Magnet[ii] = (Nodes[ii].F - Fmin) / (range);
                                                Mass[ii] = Magnet[ii] * intensity + alpha;
                                            }
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                Forces[ii].x = new double[N];
                                            }
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                int[] Neigbour = NFC.Structures(ii, NPop, structure);
                                                for (int jj = 0; jj < Neigbour.Length; jj++)
                                                {
                                                    if (ii != Neigbour[jj])
                                                    {
                                                        if (D.Text == "1")
                                                            distances[ii] = distance1(Nodes[ii].x, Nodes[Neigbour[jj]].x, N);
                                                        //else
                                                        //    distances[ii] = distance4(LHC[ii].LH, LHC[Neigbour[jj]].LH, N, P);
                                                        if (distances[ii] > short_range)
                                                        {
                                                            for (int n = 0; n < N; n++)
                                                                Forces[ii].x[n] = Forces[ii].x[n] + (Nodes[Neigbour[jj]].x[n] - Nodes[ii].x[n]) * Magnet[Neigbour[jj]] / distances[ii];

                                                        }

                                                        //else
                                                        //{
                                                        //    for (int n = 0; n < N; n++)
                                                        //            Forces[ii].x[n] = rr.NextDouble();

                                                        //    Nodes[ii].x = NFC.node_constructing(N, P, rr);
                                                        //}
                                                    }
                                                }
                                                //for (int n = 0; n < N; n++)
                                                //    for (int p = 0; p < P; p++)
                                                //    {
                                                //        Forces[ii].x[n, p] = Forces[ii].x[n, p] / (NPop - 1);
                                                //        //Forces[ii].x[n, p] = Forces[ii].x[n, p] * rr.NextDouble();
                                                //    }

                                            }
                                            if (acceleration == 0)
                                            {
                                                for (int ii = 0; ii < NPop; ii++)
                                                {
                                                    for (int n = 0; n < N; n++)
                                                    {
                                                        v[ii].x[n] = (Forces[ii].x[n] / Mass[ii]);
                                                        Nodes[ii].x[n] = Nodes[ii].x[n] + v[ii].x[n];
                                                    }

                                                    //Nodes[ii].x = NFC.Node_treatment(Nodes[ii].x, N, Rand);
                                                    Solution_Treatment(boundary_method, ref Nodes[ii].x, Rand);
                                                }
                                            }
                                            else if (acceleration == 1)
                                            {
                                                for (int ii = 0; ii < NPop; ii++)
                                                {
                                                    for (int n = 0; n < N; n++)
                                                    {
                                                        a[ii].x[n] = (Forces[ii].x[n] / Mass[ii]);
                                                        v[ii].x[n] = v[ii].x[n] * Rand.NextDouble() + a[ii].x[n];
                                                        Nodes[ii].x[n] = Nodes[ii].x[n] + v[ii].x[n];
                                                    }

                                                    Nodes[ii].x = NFC.Node_treatment(Nodes[ii].x, N, Rand);
                                                }
                                            }


                                        }
                                        costs[tt] = Best_indi.F;
                                        Nodes_test[tt].x = (double[])Best_indi.x.Clone();
                                    }

                                    sw = new StreamWriter(result);
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        sw.Write(costs[tt].ToString() + " ");
                                    }
                                    sw.Close();
                                }

                                //}

                            }
                            #endregion
                            #region PSO
                            else if (method == "PSO")
                            {
                                double W1 = Convert.ToDouble(Wmin.Text);
                                double W2 = Convert.ToDouble(Wmax.Text);
                                double w_step = Convert.ToDouble(W_step.Text);
                                double C11 = Convert.ToDouble(cmin.Text);
                                double C12 = Convert.ToDouble(cmax.Text);
                                double C11_step = Convert.ToDouble(c_step.Text);
                                double C21 = Convert.ToDouble(cmin21.Text);
                                double C22 = Convert.ToDouble(cmin22.Text);
                                double C22_step = Convert.ToDouble(c2_step.Text);
                                //double[] Cs = { 0.0001, 0.001, 0.1, 0.5, 1, 1.5 };
                                //double[] Ws = { 0.1, 0.5, 0.73, 1, 1.5 };
                                //double[] Cs = { 0.001, 1.5, 1.5, 1.5, 1.5, 0.001 };
                                //double[] Ws = { 1, 0.73, 0.73, 0.73, 0.73, 1 };
                                //for (double C1 = C11; C1 <= C12; C1 += C11_step)
                                //    for (double C2 = C21; C2 <= C12; C2 += C22_step)
                                //for (int WW = 0; WW < 5; WW++)
                                //    for (int CC = 0; CC < 6; CC++)
                                //    {
                                double[] C_s = { 0.5, 0.001, 0.001, 0.001, 0.001, 0.001, 1, 1, 0.001, 0.001, 0.001, 0.001, 0.001, 0, 0.001, 1, 1, 1, 1, 1, 0.001, 1 };
                                double[] W_s = { 0.5, 1, 1, 1, 1, 1, 0.5, 0.5, 1, 1, 1, 1, 1, 0, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1 };
                                //double C1 = Cs[OF - 22];
                                //double C2 = Cs[OF - 22];
                                //double W = Ws[OF - 22];
                                double C1 = C_s[OF];
                                double C2 = C_s[OF];
                                double W = W_s[OF];
                                //ss2 = "Results\\fit_N_" + N.ToString() + "_bound_" + boundary_method.ToString() + "_Pop_" + NPop.ToString() + "_W_" + W.ToString() + "_C1_" + C1.ToString() + "_C2_" + C2.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "OF_" + OF.ToString() + ".txt";
                                if (!File.Exists(result))
                                {
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        NF[] Nodes = new NF[NPop];
                                        NF[] BNodes = new NF[NPop];
                                        NF[] Forces = new NF[NPop];
                                        NF[] a = new NF[NPop];
                                        NF[] v = new NF[NPop];
                                        NF Best_nodes = new NF();
                                        double[] distances = new double[NPop];
                                        Best_nodes.F = -1E290;

                                        for (int ii = 0; ii < NPop; ii++)
                                        {
                                            v[ii].x = new double[N];
                                        }
                                        for (int ii = 0; ii < NPop; ii++)
                                        {
                                            Nodes[ii].x = NFC.node_constructing(N, Rand);
                                            BNodes[ii].x = (double[])Nodes[ii].x.Clone();
                                            BNodes[ii].F = Nodes[ii].F;
                                        }
                                        for (int it = 0; it < Iteration; it++)
                                        {
                                            Fmax = -1E290; int max_index = 0;
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                if (OF < 16 || OF == 22)
                                                    Nodes[ii].F = NFC.compute_fitness(Nodes[ii].x, OF, N);
                                                else
                                                    Nodes[ii].F = NFC.compute_fitness(Nodes[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                if (Fmax < Nodes[ii].F)
                                                {
                                                    Fmax = Nodes[ii].F;
                                                    max_index = ii;
                                                }
                                            }
                                            if (Best_nodes.F < Fmax)
                                            {
                                                Best_nodes.F = Fmax;
                                                Best_nodes.x = (double[])Nodes[max_index].x.Clone();

                                            }
                                            BNodes = NFC.define_pbest(Nodes, BNodes, NPop);
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                for (int n = 0; n < N; n++)
                                                {
                                                    v[ii].x[n] = W * v[ii].x[n] + C1 * Rand.NextDouble() * (BNodes[ii].x[n] - Nodes[ii].x[n]) + C2 * Rand.NextDouble() * (Best_nodes.x[n] - Nodes[ii].x[n]);
                                                    Nodes[ii].x[n] = Nodes[ii].x[n] + v[ii].x[n];
                                                }

                                                //NFC.Node_treatment(ref Nodes[ii].x, ref v[ii].x, N, Rand);
                                                Solution_Treatment(boundary_method, ref Nodes[ii].x, Rand);
                                            }
                                        }

                                        costs[tt] = Best_nodes.F;
                                        Nodes_test[tt].x = (double[])Best_nodes.x.Clone();
                                        Range.Text = C1.ToString() + " " + C2.ToString() + " " + W.ToString() + " " + OF.ToString() + " " + tt.ToString();
                                        Range.Refresh();
                                    }


                                    sw = new StreamWriter(result);
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        sw.Write(costs[tt].ToString() + " ");
                                    }
                                    sw.Close();
                                }
                                //}
                            }

                            #endregion
                            #region DE
                            else if (method == "DE")
                            {
                                //double C1 = Convert.ToDouble(Cm_ODE.Text);
                                //double C2 = Convert.ToDouble(Cmax_ODE.Text);
                                //double C_step = Convert.ToDouble(C_step_ODE.Text);
                                //double Fmi = Convert.ToDouble(Fmin_ODE.Text);
                                //double Fma = Convert.ToDouble(Fmax_ODE.Text);
                                //double F_step = Convert.ToDouble(F_s_ODE.Text);
                                //double[] Fs = { 0.001, 0.1, 0.5, 1, 2 };
                                //double[] C_s = { 0.001, 0.2, 0.4, 0.8, 1 };
                                //double[] Fs = { 0.1, 0.5, 0.1, 0.5, 0.5, 0.1 };
                                //double[] C_s = { 0.8, 0.2, 0.4, 0.2, 0.2, 0.4 };
                                //for (int II = 0; II < C_s.Length; II++)
                                //    for (int JJ = 0; JJ < Fs.Length; JJ++)
                                //    {
                                double[] Fs = { 0.1, 0.5, 0.5, 0.1, 0.1, 0.1, 0.5, 0.5, 0.1, 0.5, 0.1, 0.1, 0.1, 0, 0.1, 0.1, 0.1, 0.5, 2.0, 0.1, 0.1, 0.1 };
                                double[] C_s = { 0.2, 0.2, 0.2, 0.4, 0.2, 0.4, 0.2, 0.2, 0.4, 0.2, 0.2, 0.2, 0.2, 0, 0.4, 1, 1, 0.2, 0.2, 1, 0.4, 0.8 };
                                //double C = C_s[OF - 22];
                                //double F = Fs[OF - 22];
                                double C = C_s[OF];
                                double F = Fs[OF];
                                //ss2 = "Results\\fit_N_" + "_bound_" + boundary_method.ToString() + N.ToString() + "_Pop_" + NPop.ToString() + "_C_" + C.ToString() + "_F_" + F.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "OF_" + OF.ToString() + ".txt";
                                if (!File.Exists(result))
                                {
                                    ArrayList SQ_source = new ArrayList();
                                    for (int ii = 0; ii < NPop; ii++)
                                    {
                                        SQ_source.Add(ii);
                                    }
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        NF[] Nodes = new NF[NPop];
                                        NF[] NodesPrim = new NF[NPop];
                                        NF[] V = new NF[NPop];
                                        NF[] U = new NF[NPop];
                                        NF Best_nodes = new NF();
                                        double[] distances = new double[NPop];
                                        Best_nodes.F = -1E250;
                                        for (int ii = 0; ii < NPop; ii++)
                                        {
                                            U[ii].x = new double[N];
                                            Nodes[ii].x = NFC.node_constructing(N, Rand);
                                            if (OF < 16 || OF == 22)
                                                Nodes[ii].F = NFC.compute_fitness(Nodes[ii].x, OF, N);
                                            else
                                                Nodes[ii].F = NFC.compute_fitness(Nodes[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                        }
                                        for (int it = 0; it < Iteration; it++)
                                        {
                                            Fmax = -1E250; Fmin = 1000000000000000;
                                            int min_index = 0; int max_index = 0;
                                            for (int ii = 0; ii < NPop; ii++)
                                            {

                                                //Nodes[ii].x = NFC.Node_treatment(Nodes[ii].x, N, Rand);
                                                Solution_Treatment(boundary_method, ref Nodes[ii].x, Rand);
                                                ArrayList SQ = (ArrayList)SQ_source.Clone();
                                                V[ii].x = NFC.DE_OP(ii, Nodes, Rand, N, SQ, F);
                                                NFC.Node_treatment(V[ii].x, N, Rand);
                                                for (int jj = 0; jj < N; jj++)
                                                {
                                                    if (Rand.NextDouble() < C)
                                                        U[ii].x[jj] = V[ii].x[jj];
                                                    else
                                                        U[ii].x[jj] = Nodes[ii].x[jj];
                                                }
                                                if (OF < 16 || OF == 22)
                                                    U[ii].F = NFC.compute_fitness(U[ii].x, OF, N);
                                                else
                                                    U[ii].F = NFC.compute_fitness(U[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                if (U[ii].F >= Nodes[ii].F)
                                                {
                                                    NodesPrim[ii].x = (double[])U[ii].x.Clone();
                                                    NodesPrim[ii].F = U[ii].F;
                                                }
                                                else
                                                {
                                                    NodesPrim[ii].x = (double[])Nodes[ii].x.Clone();
                                                    NodesPrim[ii].F = Nodes[ii].F;
                                                }
                                                if (Fmax < Nodes[ii].F)
                                                {
                                                    Fmax = Nodes[ii].F;
                                                    max_index = ii;
                                                }
                                            }
                                            if (Best_nodes.F < Fmax)
                                            {
                                                Best_nodes.F = Fmax;
                                                Best_nodes.x = (double[])Nodes[min_index].x.Clone();

                                            }
                                            NFC.PprimtoP(ref Nodes, ref NodesPrim);

                                        }

                                        costs[tt] = Best_nodes.F;
                                        Nodes_test[tt].x = (double[])Best_nodes.x.Clone();
                                        Range.Text = C.ToString() + " " + F.ToString() + " " + OF.ToString() + " " + tt.ToString();
                                        Range.Refresh();
                                    }
                                    sw = new StreamWriter(result);
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        sw.Write(costs[tt].ToString() + " ");
                                    }
                                    sw.Close();
                                }
                                //}
                            }

                            #endregion
                            #region ODE
                            else if (method == "ODE")
                            {
                                //double[] Fs = { 0.2000, 0.2000, 0.2000, 0.2000, 0.2000, 0.2000, 0.2000, 0.2000, 0.2000, 0.2000, 0.2000, 0.2000, 0.2000, 0.2 };
                                //double[] C_s = { 0.4, 0.4, 0.4, 0.4, 0.4, 0.4 };
                                //double[] F_s = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 };
                                double Jr = 0.3;
                                //for (int F_index = 0; F_index < 5; F_index++)
                                //    for (int C_index = 0; C_index < 5; C_index++)
                                //    {

                                //double C = C_s[OF - 22];
                                //double F = F_s[OF - 22];
                                double C = 0.4;
                                double F = 0.1;
                                //ss2 = "Results\\fit_N_" + "_bound_" + boundary_method.ToString() + N.ToString() + "_Pop_" + NPop.ToString() + "_C_" + C.ToString() + "_F_" + F.ToString() + "_Jr_" + Jr.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "OF_" + OF.ToString() + ".txt";
                                if (!File.Exists(result))
                                {
                                    ArrayList SQ_source = new ArrayList();
                                    for (int ii = 0; ii < NPop; ii++)
                                    {
                                        SQ_source.Add(ii);
                                    }
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        NF[] Nodes = new NF[NPop];
                                        NF[] NodesPrim = new NF[NPop];
                                        NF[] ONodes = new NF[NPop];
                                        NF[] V = new NF[NPop];
                                        NF[] U = new NF[NPop];
                                        NF Best_nodes = new NF();
                                        double[] distances = new double[NPop];
                                        Best_nodes.F = -1E250;
                                        for (int ii = 0; ii < NPop; ii++)
                                        {
                                            U[ii].x = new double[N];
                                            Nodes[ii].x = NFC.node_constructing(N, Rand);
                                            ONodes[ii].x = NFC.Opposition_Learning(Nodes[ii].x);
                                            double temp1 = 0; double temp2 = 0;
                                            if (OF < 16 || OF == 22)
                                                temp1 = NFC.compute_fitness(Nodes[ii].x, OF, N);
                                            else
                                                temp1 = NFC.compute_fitness(Nodes[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                            if (OF < 16 || OF == 22)
                                                temp2 = NFC.compute_fitness(ONodes[ii].x, OF, N);
                                            else
                                                temp2 = NFC.compute_fitness(ONodes[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                            if (temp1 >= temp2)
                                                Nodes[ii].F = temp1;
                                            else
                                            {
                                                Nodes[ii].x = (double[])ONodes[ii].x.Clone();
                                                Nodes[ii].F = temp2;
                                            }
                                        }
                                        for (int it = 0; it < Iteration; it++)
                                        {
                                            Fmax = -1E250; Fmin = 1E250;
                                            int min_index = 0; int max_index = 0;
                                            for (int ii = 0; ii < NPop; ii++)
                                            {

                                                //Nodes[ii].x = NFC.Node_treatment(Nodes[ii].x, N, Rand);
                                                Solution_Treatment(boundary_method, ref Nodes[ii].x, Rand);
                                                ArrayList SQ = (ArrayList)SQ_source.Clone();
                                                V[ii].x = NFC.DE_OP(ii, Nodes, Rand, N, SQ, F);
                                                if (OF < 16 || OF == 22)
                                                    Nodes[ii].F = NFC.compute_fitness(Nodes[ii].x, OF, N);
                                                else
                                                    Nodes[ii].F = NFC.compute_fitness(Nodes[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                for (int jj = 0; jj < N; jj++)
                                                {
                                                    if (Rand.NextDouble() < C)
                                                        U[ii].x[jj] = V[ii].x[jj];
                                                    else
                                                        U[ii].x[jj] = Nodes[ii].x[jj];
                                                }
                                                if (OF < 16 || OF == 22)
                                                    U[ii].F = NFC.compute_fitness(U[ii].x, OF, N);
                                                else
                                                    U[ii].F = NFC.compute_fitness(U[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                if (U[ii].F >= Nodes[ii].F)
                                                {
                                                    NodesPrim[ii].x = (double[])U[ii].x.Clone();
                                                    NodesPrim[ii].F = U[ii].F;
                                                }
                                                else
                                                {
                                                    NodesPrim[ii].x = (double[])Nodes[ii].x.Clone();
                                                    NodesPrim[ii].F = Nodes[ii].F;
                                                }
                                                if (Fmax < Nodes[ii].F)
                                                {
                                                    Fmax = Nodes[ii].F;
                                                    max_index = ii;
                                                }
                                                if (Fmin > Nodes[ii].F)
                                                {
                                                    Fmin = Nodes[ii].F;
                                                    min_index = ii;
                                                }
                                            }
                                            if (Best_nodes.F < Fmax)
                                            {
                                                Best_nodes.F = Fmax;
                                                Best_nodes.x = (double[])Nodes[min_index].x.Clone();

                                            }
                                            NFC.PprimtoP(ref Nodes, ref NodesPrim);
                                            if (Rand.NextDouble() < Jr)
                                            {
                                                for (int ii = 0; ii < NPop; ii++)
                                                {
                                                    ONodes[ii].x = NFC.generating_jumping_orignal(Nodes, Nodes[ii].x, NPop, N, ii);
                                                    if (OF < 16 || OF == 22)
                                                        ONodes[ii].F = NFC.compute_fitness(ONodes[ii].x, OF, N);
                                                    else
                                                        ONodes[ii].F = NFC.compute_fitness(ONodes[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                }
                                                Nodes = NFC.define_pbest(Nodes, ONodes, NPop);
                                            }

                                        }

                                        costs[tt] = Best_nodes.F;
                                        Nodes_test[tt].x = (double[])Best_nodes.x.Clone();
                                        Range.Text = C.ToString() + " " + F.ToString() + " " + Jr.ToString() + " " + OF.ToString() + " " + tt.ToString();
                                        Range.Refresh();
                                    }

                                    sw = new StreamWriter(result);
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        sw.Write(costs[tt].ToString() + " ");
                                    }
                                    sw.Close();
                                }
                                //}
                            }

                            #endregion
                            #region ES
                            else if (method == "ES")
                            {
                                //double[] Lambda_s = { 0.2, 0.4, 0.6, 0.8, 1 };
                                ////// sigma constant which is the constant factor to be multiplied to sigma
                                //double[] SC_s = { 1.01, 1.1, 1.2, 1.5, 2 };
                                //double[] Lambda_s = { 1, 0.8, 0.6, 1, 1, 0.8 };
                                //// sigma constant which is the constant factor to be multiplied to sigma
                                //double[] SC_s = { 1.01, 1.01, 1.01, 1.01, 1.01, 1.01 };
                                //for (int cc = 0; cc < 5; cc++)
                                //    for (int ww = 0; ww < 5; ww++)
                                //    {
                                double[] Lambda_s = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.6 };
                                double[] SC_s = { 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.5, 1.1, 1.1, 1.2, 5, 1.1, 1.1 };
                                int breeds_number = Convert.ToInt32(Lambda_s[OF ] * NPop);
                                double SC = SC_s[OF];
                                //ss = "Result_N_" + N.ToString() + "_Pop_" + NPop.ToString() + "_Lambda_" + Lambda_s[cc].ToString() + "_Sigma_Fac_" + SC.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_OF" + OF.ToString() + ".txt";
                                //ss2 = "Results\\fit_N_" + "_bound_" + boundary_method.ToString() + N.ToString() + "_Pop_" + NPop.ToString() + "_Lambda_" + Lambda_s[OF].ToString() + "_Sigma_Fac_" + SC.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_OF" + OF.ToString() + ".txt";
                                if (!File.Exists(result))
                                {
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        NF[] Pop = new NF[NPop];
                                        NF Best_nodes = new NF();
                                        Best_nodes.F = -1E250;
                                        NF[] Breeds = new NF[breeds_number];
                                        double[] distances = new double[NPop];
                                        int phi_count = 0;
                                        double sigma = 1;
                                        for (int ii = 0; ii < NPop; ii++)
                                        {
                                            Pop[ii].x = NFC.node_constructing(N, Rand);
                                            if (OF < 16 || OF == 22)
                                                Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N);
                                            else
                                                Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                        }

                                        for (int it = 0; it < Iteration; it++)
                                        {
                                            NF[] All = new NF[NPop + breeds_number];
                                            for (int ii = 0; ii < breeds_number; ii++)
                                            {
                                                Breeds[ii].x = new double[N];
                                                int parent1 = Rand.Next(NPop);
                                                int parent2 = Rand.Next(NPop);
                                                for (int n = 0; n < N; n++)
                                                {
                                                    double rand = Rand.NextDouble();
                                                    if (rand < 0.5)
                                                    {
                                                        Breeds[ii].x[n] = Pop[parent1].x[n];
                                                    }
                                                    else
                                                    {
                                                        Breeds[ii].x[n] = Pop[parent2].x[n];
                                                    }
                                                }
                                            }
                                            for (int ii = 0; ii < breeds_number; ii++)
                                            {
                                                for (int n = 0; n < N; n++)
                                                {
                                                    Breeds[ii].x[n] = Breeds[ii].x[n] + sigma * NFC.Normal_Distribution(Rand, 0, 1);
                                                }
                                                //NFC.Node_treatment(Breeds[ii].x, N, Rand);
                                                Solution_Treatment(boundary_method, ref Breeds[ii].x, Rand);
                                                if (OF < 16 || OF == 22)
                                                    Breeds[ii].F = NFC.compute_fitness(Breeds[ii].x, OF, N);
                                                else
                                                    Breeds[ii].F = NFC.compute_fitness(Breeds[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                All[ii].x = (double[])Breeds[ii].x.Clone();
                                                All[ii].F = Breeds[ii].F;

                                            }

                                            double[] OldCost = new double[NPop];
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                if (OF < 16 || OF == 22)
                                                    OldCost[ii] = NFC.compute_fitness(Pop[ii].x, OF, N);
                                                else
                                                    OldCost[ii] = NFC.compute_fitness(Pop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                for (int n = 0; n < N; n++)
                                                {
                                                    Pop[ii].x[n] = Pop[ii].x[n] + sigma * NFC.Normal_Distribution(Rand, 0, 1);
                                                }
                                                NFC.Node_treatment(Pop[ii].x, N, Rand);
                                                if (OF < 16 || OF == 22)
                                                    Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N);
                                                else
                                                    Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);

                                                // putting all pop into a temporary combinator
                                                All[ii + breeds_number].x = (double[])Pop[ii].x.Clone();
                                                All[ii + breeds_number].F = Pop[ii].F;

                                                if (Pop[ii].F > OldCost[ii])
                                                    phi_count++;
                                            }
                                            int mutation_count = it * NPop;
                                            if (phi_count > (mutation_count / 5))
                                                sigma = sigma / SC;
                                            else if (phi_count < (mutation_count / 5))
                                                sigma = sigma * SC;
                                            // Sort all POP+Breed
                                            All = NFC.mergesort(All, All.Length);

                                            //storing the best node into best node container
                                            if (Best_nodes.F < All[0].F)
                                            {
                                                Best_nodes.x = (double[])All[0].x.Clone();
                                                Best_nodes.F = All[0].F;
                                            }

                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                Pop[ii].x = (double[])All[ii].x.Clone();
                                                Pop[ii].F = All[ii].F;
                                            }
                                        }

                                        costs[tt] = Best_nodes.F;
                                        Nodes_test[tt].x = (double[])Best_nodes.x.Clone();

                                    }

                                    //sw = new StreamWriter(ss);
                                    //for (int tt = 0; tt < Ntest; tt++)
                                    //{
                                    //    for (int ii = 0; ii < N; ii++)
                                    //    {
                                    //            sw.Write(Nodes_test[tt].x[ii].ToString() + " ");
                                    //    }
                                    ////}
                                    //sw.Close();
                                    sw = new StreamWriter(result);
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        sw.Write(costs[tt].ToString() + " ");
                                    }
                                    sw.Close();
                                }

                                //}

                            }

                            #endregion
                            #region FES
                            else if (method == "FES")
                            {

                                //double[] Lambda_s = { 0.2, 0.4, 0.6, 0.8, 1 };
                                ////// sigma constant which is the constant factor to be multiplied to sigma
                                //double[] SC_s = { 1.1, 1.2, 1.5, 2, 5 };
                                //double[] Lambda_s = { 1, 1, 1, 1, 1, 1 };
                                //// sigma constant which is the constant factor to be multiplied to sigma
                                //double[] SC_s = { 1.1, 1.1, 1.1, 1.1, 1.1, 1.1 };
                                //for (int cc = 0; cc < 5; cc++)
                                //    for (int ww = 0; ww < 5; ww++)
                                //    {
                                double[] Lambda_s = { 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 0, 1.0000, 0.8000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 };
                                double[] SC_s = { 1.1000, 1.1000, 1.1000, 1.1000, 1.1000, 1.1000, 1.1000, 1.1000, 1.1000, 1.1000, 1.1000, 1.1000, 1.1000, 0, 1.1000, 2.0000, 1.1000, 1.1000, 1.2000, 1.1000, 1.1000, 1.1000 };
                                int breeds_number = Convert.ToInt32(Lambda_s[OF ] * NPop);
                                double SC = SC_s[OF];
                                //ss = "Result_N_" + N.ToString() + "_Pop_" + NPop.ToString() + "_Lambda_" + Lambda_s[cc].ToString() + "_Sigma_Fac_" + SC.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_OF" + OF.ToString() + ".txt";
                                //ss2 = "Results\\fit_N_" + N.ToString() + "_bound_" + boundary_method.ToString() + "_Pop_" + NPop.ToString() + "_Lambda_" + Lambda_s[OF].ToString() + "_Sigma_Fac_" + SC.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_OF" + OF.ToString() + ".txt";
                                if (!File.Exists(result))
                                {
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        NF[] Pop = new NF[NPop];
                                        NF Best_nodes = new NF();
                                        Best_nodes.F = -1E250;
                                        NF[] Breeds = new NF[breeds_number];
                                        double[] distances = new double[NPop];
                                        int phi_count = 0;
                                        double sigma = 1;
                                        for (int ii = 0; ii < NPop; ii++)
                                        {
                                            Pop[ii].x = NFC.node_constructing(N, Rand);
                                            if (OF < 16 || OF == 22)
                                                Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N);
                                            else
                                                Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                        }

                                        for (int it = 0; it < Iteration; it++)
                                        {
                                            NF[] All = new NF[NPop + breeds_number];
                                            for (int ii = 0; ii < breeds_number; ii++)
                                            {
                                                Breeds[ii].x = new double[N];
                                                int parent1 = Rand.Next(NPop);
                                                int parent2 = Rand.Next(NPop);
                                                for (int n = 0; n < N; n++)
                                                {
                                                    double rand = Rand.NextDouble();
                                                    if (rand < 0.5)
                                                    {
                                                        Breeds[ii].x[n] = Pop[parent1].x[n];
                                                    }
                                                    else
                                                    {
                                                        Breeds[ii].x[n] = Pop[parent2].x[n];
                                                    }
                                                }
                                            }
                                            for (int ii = 0; ii < breeds_number; ii++)
                                            {
                                                for (int n = 0; n < N; n++)
                                                {
                                                    Breeds[ii].x[n] = Breeds[ii].x[n] + sigma * NFC.Cauchy(Rand, 0, 1);
                                                }
                                                //NFC.Node_treatment(Breeds[ii].x, N, Rand);
                                                Solution_Treatment(boundary_method,ref Breeds[ii].x, Rand);
                                                if (OF < 16 || OF == 22)
                                                    Breeds[ii].F = NFC.compute_fitness(Breeds[ii].x, OF, N);
                                                else
                                                    Breeds[ii].F = NFC.compute_fitness(Breeds[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                All[ii].x = (double[])Breeds[ii].x.Clone();
                                                All[ii].F = Breeds[ii].F;

                                            }

                                            double[] OldCost = new double[NPop];
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                if (OF < 16 || OF == 22)
                                                    OldCost[ii] = NFC.compute_fitness(Pop[ii].x, OF, N);
                                                else
                                                    OldCost[ii] = NFC.compute_fitness(Pop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                for (int n = 0; n < N; n++)
                                                {
                                                    Pop[ii].x[n] = Pop[ii].x[n] + sigma * NFC.Cauchy(Rand, 0, 1);
                                                }
                                                NFC.Node_treatment(Pop[ii].x, N, Rand);
                                                if (OF < 16 || OF == 22)
                                                    Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N);
                                                else
                                                    Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);

                                                // putting all pop into a temporary combinator
                                                All[ii + breeds_number].x = (double[])Pop[ii].x.Clone();
                                                All[ii + breeds_number].F = Pop[ii].F;

                                                if (Pop[ii].F > OldCost[ii])
                                                    phi_count++;
                                            }
                                            int mutation_count = it * NPop;
                                            if (phi_count > (mutation_count / 5))
                                                sigma = sigma / SC;
                                            else if (phi_count < (mutation_count / 5))
                                                sigma = sigma * SC;
                                            // Sort all POP+Breed
                                            All = NFC.mergesort(All, All.Length);

                                            //storing the best node into best node container
                                            if (Best_nodes.F < All[0].F)
                                            {
                                                Best_nodes.x = (double[])All[0].x.Clone();
                                                Best_nodes.F = All[0].F;
                                            }

                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                Pop[ii].x = (double[])All[ii].x.Clone();
                                                Pop[ii].F = All[ii].F;
                                            }
                                        }

                                        costs[tt] = Best_nodes.F;
                                        Nodes_test[tt].x = (double[])Best_nodes.x.Clone();

                                    }
                                    sw = new StreamWriter(result);
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        sw.Write(costs[tt].ToString() + " ");
                                    }
                                    sw.Close();
                                }
                                //}

                            }

                            #endregion
                            #region EP
                            else if (method == "EP")
                            {
                                //double[] Tournament_s = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1 };
                                //double[] Tournament_s = { 1, 0.5, 0.1, 0.5, 0.3, 0.5 };
                                //double[] Tournament_s = new double[25];
                                //double dist = 1.0 / 25;
                                //Tournament_s[0] = dist;
                                //for (int ii = 1; ii < 25; ii++)
                                //{
                                //    Tournament_s[ii] = dist + Tournament_s[ii - 1];
                                //}
                                double[] etha = new double[N];
                                double[] Tournament_s = { 0.3, .2, 0.2, 0.4, .2, .3, 0.6, 0.1, 0.1, 0.1, 0.1, 0.2, 0.1, 0, 0.1, 0.2, .1, 0.1, 0.9, 0.1, 0.3, 0.7 };
                                double[] etha_prim = new double[N];
                                // sigma constant which is the constant factor to be multiplied to sigma
                                //double[] SC_s = { 1.1, 1.2, 1.5, 2, 5 };
                                //for (int cc = 0; cc < Tournament_s.Length; cc++)
                                //{
                                int tournament_size = Convert.ToInt32(Math.Ceiling(Tournament_s[OF] * NPop));
                                //ss2 = "Results\\fit_N_" + N.ToString() + "_bound_" + boundary_method.ToString() + "_Pop_" + NPop.ToString() + "_Tournam_" + tournament_size.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_OF" + OF.ToString() + ".txt";
                                if (!File.Exists(result))
                                {
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        NF[] Pop = new NF[NPop];
                                        NF Best_nodes = new NF();
                                        Best_nodes.F = -1E250;
                                        NF[] Breeds = new NF[NPop];
                                        NF[] All = new NF[NPop * 2];
                                        double sigma = 1 / (Math.Sqrt(2.0 * Math.Sqrt(N)));
                                        double sigma_prim = 1 / (Math.Sqrt(2.0 * N));
                                        for (int j = 0; j < N; j++)
                                        {
                                            etha[j] = Rand.NextDouble(); etha_prim[j] = Rand.NextDouble();
                                        }
                                        // Initialization
                                        for (int ii = 0; ii < NPop; ii++)
                                        {
                                            Pop[ii].x = NFC.node_constructing(N, Rand);

                                            All[ii].x = (double[])Pop[ii].x.Clone();
                                            if (OF < 16 || OF == 22)
                                                Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N);
                                            else
                                                Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                            All[ii].F = Pop[ii].F;
                                            if (Best_nodes.F < Pop[ii].F)
                                            {
                                                Best_nodes.F = Pop[ii].F;
                                                Best_nodes.x = (double[])Pop[ii].x.Clone();
                                            }
                                        }
                                        ///////////////////////////////

                                        for (int it = 0; it < Iteration; it++)
                                        {
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                Breeds[ii].x = new double[N];
                                                for (int jj = 0; jj < N; jj++)
                                                {
                                                    Breeds[ii].x[jj] = Pop[ii].x[jj] + etha[jj] * NFC.Normal_Distribution(Rand, 0, 1);
                                                    etha_prim[jj] = etha[jj] * (Math.Exp(sigma_prim * NFC.Normal_Distribution(Rand, 0, 1) + sigma * NFC.Normal_Distribution(Rand, 0, 1)));
                                                }
                                                //NFC.Node_treatment(Breeds[ii].x, N, Rand);
                                                Solution_Treatment(boundary_method, ref Breeds[ii].x, Rand);
                                                All[NPop + ii].x = (double[])Breeds[ii].x.Clone();
                                                if (OF < 16 || OF == 22)
                                                    Breeds[ii].F = NFC.compute_fitness(Breeds[ii].x, OF, N);
                                                else
                                                    Breeds[ii].F = NFC.compute_fitness(Breeds[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                All[NPop + ii].F = Breeds[ii].F;
                                                if (Best_nodes.F < Breeds[ii].F)
                                                {
                                                    Best_nodes.F = Breeds[ii].F;
                                                    Best_nodes.x = (double[])Breeds[ii].x;
                                                }
                                            }

                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                int[] Tournament_indi = NFC.sq_rand_gen2(tournament_size, NPop * 2, Rand);
                                                double max = -1E300; int max_index = -1;
                                                for (int jj = 0; jj < tournament_size; jj++)
                                                {
                                                    if (max < All[Tournament_indi[jj]].F)
                                                    {
                                                        max = All[Tournament_indi[jj]].F;
                                                        max_index = Tournament_indi[jj];
                                                    }
                                                }
                                                Pop[ii].x = (double[])All[max_index].x.Clone();
                                                Pop[ii].F = max;
                                            }
                                            // transferring the survivor to their places
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                All[ii].x = (double[])Pop[ii].x.Clone();
                                                All[ii].F = Pop[ii].F;
                                            }
                                        }

                                        costs[tt] = Best_nodes.F;
                                        Nodes_test[tt].x = (double[])Best_nodes.x.Clone();

                                    }
                                    sw = new StreamWriter(result);
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        sw.Write(costs[tt].ToString() + " ");
                                    }
                                    sw.Close();
                                }
                                //}

                            }

                            #endregion
                            #region FEP
                            else if (method == "FEP")
                            {
                                //double[] Tournament_s = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1 };
                                //double[] Tournament_s = { 0.4, 0.5, 0.8, 0.2, 1, 0.4 };
                                //double[] Tournament_s = new double[25];
                                //double dist = 1.0 / 25;
                                //Tournament_s[0] = dist;
                                //for (int ii = 1; ii < 25; ii++)
                                //{
                                //    Tournament_s[ii] = dist + Tournament_s[ii - 1];
                                //}
                                double[] Tournament_s = { 1.0, 0.2, .5, 0.9, 0.1, 0.4, 0.4, 0.4, 0.7, 0.7, 0.1, 0.5, 0.4, 0, 0.6, 0.1, 1.0, 0.1, 0.1, 0.5, 0.6, 1 };
                                double[] etha = new double[N];
                                double[] etha_prim = new double[N];
                                //for (int cc = 0; cc < Tournament_s.Length; cc++)
                                //{
                                int tournament_size = Convert.ToInt32(Math.Ceiling(Tournament_s[OF] * NPop));
                                //ss2 = "Results\\fit_N_" + N.ToString() + "_bound_" + boundary_method.ToString() + "_Pop_" + NPop.ToString() + "_Tournam_" + tournament_size.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_OF" + OF.ToString() + ".txt";
                                if (!File.Exists(result))
                                {
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        NF[] Pop = new NF[NPop];
                                        NF Best_nodes = new NF();
                                        Best_nodes.F = -1E280;
                                        NF[] Breeds = new NF[NPop];
                                        NF[] All = new NF[NPop * 2];
                                        double sigma = 1 / (Math.Sqrt(2.0 * Math.Sqrt(N)));
                                        double sigma_prim = 1 / (Math.Sqrt(2.0 * N));
                                        for (int j = 0; j < N; j++)
                                        {
                                            etha[j] = Rand.NextDouble(); etha_prim[j] = Rand.NextDouble();
                                        }
                                        // Initialization
                                        for (int ii = 0; ii < NPop; ii++)
                                        {
                                            Pop[ii].x = NFC.node_constructing(N, Rand);

                                            All[ii].x = (double[])Pop[ii].x.Clone();
                                            if (OF < 16 || OF == 22)
                                                Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N);
                                            else
                                                Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                            All[ii].F = Pop[ii].F;
                                            if (Best_nodes.F < Pop[ii].F)
                                            {
                                                Best_nodes.F = Pop[ii].F;
                                                Best_nodes.x = (double[])Pop[ii].x.Clone();
                                            }
                                        }
                                        ///////////////////////////////

                                        for (int it = 0; it < Iteration; it++)
                                        {
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                Breeds[ii].x = new double[N];
                                                for (int jj = 0; jj < N; jj++)
                                                {
                                                    Breeds[ii].x[jj] = Pop[ii].x[jj] + etha[jj] * NFC.Cauchy(Rand, 0, 1);
                                                    etha_prim[jj] = etha[jj] * (Math.Exp(sigma_prim * NFC.Normal_Distribution(Rand, 0, 1) + sigma * NFC.Normal_Distribution(Rand, 0, 1)));
                                                }
                                                //NFC.Node_treatment(Breeds[ii].x, N, Rand);
                                                Solution_Treatment(boundary_method, ref Breeds[ii].x, Rand);
                                                All[NPop + ii].x = (double[])Breeds[ii].x.Clone();
                                                if (OF < 16 || OF == 22)
                                                    Breeds[ii].F = NFC.compute_fitness(Breeds[ii].x, OF, N);
                                                else
                                                    Breeds[ii].F = NFC.compute_fitness(Breeds[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                All[NPop + ii].F = Breeds[ii].F;
                                                if (Best_nodes.F < Breeds[ii].F)
                                                {
                                                    Best_nodes.F = Breeds[ii].F;
                                                    Best_nodes.x = (double[])Breeds[ii].x;
                                                }
                                            }

                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                int[] Tournament_indi = NFC.sq_rand_gen2(tournament_size, NPop * 2, Rand);
                                                double max = -1E299; int max_index = -1;
                                                for (int jj = 0; jj < tournament_size; jj++)
                                                {
                                                    if (max < All[Tournament_indi[jj]].F)
                                                    {
                                                        max = All[Tournament_indi[jj]].F;
                                                        max_index = Tournament_indi[jj];
                                                    }
                                                }
                                                Pop[ii].x = (double[])All[max_index].x.Clone();
                                                Pop[ii].F = max;
                                            }
                                            // transferring the survivor to their places
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                All[ii].x = (double[])Pop[ii].x.Clone();
                                                All[ii].F = Pop[ii].F;
                                            }
                                        }

                                        costs[tt] = Best_nodes.F;
                                        Nodes_test[tt].x = (double[])Best_nodes.x.Clone();

                                    }
                                    sw = new StreamWriter(result);
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        sw.Write(costs[tt].ToString() + " ");
                                    }
                                    sw.Close();
                                }
                                //}

                            }
                            #endregion
                            #region MA-SW
                            else if (method == "MA-SW")
                            {
                               O = NPop * Convert.ToInt32(NI.Text);
                                //double[] N_LS_GL = { 0.9, 0.1, 0.5, 0.1, 0.1, 0.9 };
                                //double[] Mu_rate_s = { 0.01, 0.01, 0.05, 0.001, 0.001, 0.01 };
                                //for (int Portion_index = 0; Portion_index < 5; Portion_index++)
                                //    for (int M_index = 0; M_index < 5; M_index++)
                                //    {
                                double[] N_LS_GL = { 0.3, 0.5, 0.5, 0.3, 0.5, 0.5, 0.5, 0.5, 0.3, 0.5, 0.7, 0.3, 0.5, 0, 0.3, 0.5, 0.5, 0.5, 0.5, 0.7, 0.3, 0.5 };
                                double[] Mu_rate_s = { 0.00001, 0.5, 0.001, 0.001, 0.00001, 5, 0.5, 0.5, 0.001, 0.01, 0.001, 0.01, 0.5, 0, 0.001, 0.1, 0.1, 0.5, 0.5, 0.00001, 0.00001, 0.1 };
                                double LS_LG = N_LS_GL[OF ];
                                double Mu_rate = Mu_rate_s[OF ];
                                Iteration = Convert.ToInt32((O * (1 - LS_LG)) / NPop) + 1;
                                int I_str = Convert.ToInt32(O * LS_LG / Iteration);
                                int best_index = 0;
                                int local_search = 0;
                                //if (Mu_rate == 0.00001)
                                //    ss2 = "fit_N_" + N.ToString() + "_Mr_" + Mu_rate.ToString() + "LSLG" + LS_LG.ToString() + "_Pop_" + NPop.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_OF" + OF.ToString() + ".txt";
                                //else
                                //ss2 = "Results\\fit_N_" + N.ToString() + "_bound_" + boundary_method.ToString() + "_Mr_" + Mu_rate.ToString() + "LSLG" + LS_LG.ToString() + "_Pop_" + NPop.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_OF" + OF.ToString() + ".txt";
                                if (!File.Exists(result))
                                {
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        LS[] Pop = new LS[NPop];
                                        NF Best_nodes = new NF();
                                        NF[] Breed = new NF[NPop];
                                        //NF[] Breed = new NF[2*NPop];//for laplacian operator
                                        Best_nodes.F = -1E250;
                                        //LS[] LS_S = new LS[N];
                                        for (int ii = 0; ii < NPop; ii++)
                                        {
                                            Pop[ii].x = NFC.node_constructing(N, Rand);

                                            if (OF < 16 || OF == 22)
                                                Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N);
                                            else
                                                Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                            Pop[ii].local_search = true;

                                        }
                                        Pop = NFC.mergesort(Pop, NPop);
                                        for (int it = 0; it < Iteration; it++)
                                        {
                                            //perform Steady State GA
                                            #region generate Breeds
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                if (Mu_rate > Rand.NextDouble())
                                                    NFC.Mutation_BGA(ref Pop[ii], Pop, N, Rand);
                                                if (OF < 16 || OF == 22)
                                                    Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N);
                                                else
                                                    Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                            }
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                int parent1 = Rand.Next(NPop);
                                                Breed[ii].x = new double[N];
                                                //Breed[ii+1].x = new double[N];
                                                //int[] S = NFC.Structures(ii, NPop, Struct);
                                                int parent2 = NFC.negative_assortative_mating_strategy(Pop[parent1].x, Pop, Rand);
                                                Breed[ii].x = NFC.BLX_Crossover(Pop[parent1].x, Pop[parent2].x);
                                                //laplac temp = NFC.Laplacian_recombination(Pop[parent1], Pop[parent2], N, Rand, 0, 0.5);
                                                //Breed[ii].x = (double[])temp.par1.Clone();
                                                //Breed[ii+1].x = (double[])temp.par2.Clone();
                                                //Breed[ii].x = NFC.Node_treatment(Breed[ii].x, N, Rand);
                                                Solution_Treatment(boundary_method, ref Breed[ii].x, Rand);
                                                //Breed[ii + 1].x = NFC.Node_treatment(Breed[ii + 1].x, N, Rand);
                                                if (OF < 16 || OF == 22)
                                                    Breed[ii].F = NFC.compute_fitness(Breed[ii].x, OF, N);
                                                else
                                                    Breed[ii].F = NFC.compute_fitness(Breed[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                //if (OF < 16 || OF == 22)
                                                //    Breed[ii+1].F = NFC.compute_fitness(Breed[ii+1].x, OF, N);
                                                //else
                                                //    Breed[ii+1].F = NFC.compute_fitness(Breed[ii+1].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                            }
                                            #endregion
                                            // Replacement Strategy

                                            Pop = NFC.mergesort(Pop, NPop);
                                            Breed = NFC.mergesort(Breed, NPop);
                                            if (Breed[0].F > Pop[NPop - 1].F)
                                            {
                                                Pop[NPop - 1].x = (double[])Breed[0].x.Clone();
                                                Pop[NPop - 1].F = Breed[0].F;
                                                Pop[NPop - 1].local_search = true;
                                            }
                                            //
                                            //Pick the best individual of Population
                                            if (local_search == NPop - 1)
                                            {
                                                for (int ii = 1; ii < NPop; ii++)
                                                {
                                                    Pop[ii].x = NFC.node_constructing(N, Rand);
                                                    if (OF < 16 || OF == 22)
                                                        Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N);
                                                    else
                                                        Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                    Pop[ii].local_search = true;
                                                    Pop[ii].visited = 0;
                                                }
                                                local_search = 0;
                                            }

                                            LS Best_Indi_BLS = new LS();
                                            Pop = NFC.mergesort(Pop, NPop);
                                            for (int ii = 0; ii < NPop; ii++)
                                                if (Pop[ii].local_search)
                                                {
                                                    Best_Indi_BLS.F = Pop[ii].F;
                                                    Best_Indi_BLS.x = (double[])Pop[ii].x.Clone();
                                                    Best_Indi_BLS.visited = Pop[ii].visited;
                                                    Best_Indi_BLS.rho = Pop[ii].rho;
                                                    if (Pop[ii].visited == 1)
                                                        Best_Indi_BLS.bias = (double[])Pop[ii].bias.Clone();
                                                    best_index = ii;
                                                    break;
                                                }
                                            // Initialize Bias and Rho
                                            if (Best_Indi_BLS.visited == 0)
                                            {
                                                Best_Indi_BLS.bias = new double[N];
                                                Best_Indi_BLS.rho = 0.01;
                                            }
                                            // Local Search Operator
                                            LS Best_Indi_ALS = new LS();
                                            Best_Indi_ALS.x = (double[])Best_Indi_BLS.x.Clone();
                                            Best_Indi_ALS.F = Best_Indi_BLS.F;
                                            Best_Indi_ALS.rho = Best_Indi_BLS.rho;
                                            Best_Indi_ALS.visited = Best_Indi_BLS.visited;
                                            Best_Indi_ALS.bias = (double[])Best_Indi_BLS.bias.Clone();
                                            NFC.SolisWets(ref Best_Indi_ALS.x, ref Best_Indi_ALS.F, N, ref Best_Indi_ALS.bias, ref Best_Indi_ALS.rho, I_str, Rand, OF, Perm, Orthogonal_Matrix, M, o);

                                            // Remove Unchanged particle from population
                                            if (Best_Indi_BLS.F - Best_Indi_ALS.F == 0)
                                            {
                                                Pop[best_index].local_search = false;
                                                //  Pop[best_index].visited = 1;
                                                local_search++;
                                            }
                                            else
                                            {
                                                Pop[best_index].x = (double[])Best_Indi_ALS.x.Clone();
                                                Pop[best_index].visited = 1;
                                                Pop[best_index].F = Best_Indi_ALS.F;
                                                Pop[best_index].rho = Best_Indi_ALS.rho;
                                                Pop[best_index].bias = (double[])Best_Indi_ALS.bias.Clone();
                                                Pop[best_index].local_search = true;
                                            }

                                            Pop = NFC.mergesort(Pop, NPop);
                                            if (Best_nodes.F < Pop[0].F)
                                            {
                                                Best_nodes.F = Pop[0].F;
                                                Best_nodes.x = (double[])Pop[0].x.Clone();
                                            }

                                            /////////////////////

                                        }
                                        costs[tt] = Best_nodes.F;

                                    }


                                    sw = new StreamWriter(result);
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        sw.Write(costs[tt].ToString() + " ");
                                    }
                                    sw.Close();
                                }
                                //}

                            }

                            #endregion
                            #region SRE
                            if (method == "SRE")
                            {
                                double[] Crossover_methods = new double[3];
                                int min_index = 0;
                                double mu1 = Convert.ToDouble(Mu1.Text);
                                double mu2 = Convert.ToDouble(Mu2.Text);
                                double mu_step = Convert.ToDouble(Mu_step.Text);
                                double[] muindex = { 0.002, 0.003, 0.004, 0.005, 0.01 };
                                double[] SP = new double[3];
                                double cro1 = Convert.ToDouble(cr1.Text);
                                double cro2 = Convert.ToDouble(cr2.Text);
                                double cro_step = Convert.ToDouble(cr_step.Text);
                                int[] Thresold_s = { 0, 1, 10, 50, 100 };
                                double[] Reg_Par = { 0.95, 0.97, 0.99 };
                                //for (double MM = mu1; MM <= mu2; MM += mu_step)
                                //    for (double cr_rate = cro1; cr_rate <= cro2; cr_rate += cro_step)
                                for (int RP = 0; RP < Reg_Par.Length; RP++)
                                    for (int Tind = 0; Tind < Thresold_s.Length; Tind++)
                                    {

                                        double cr_rate = 1;
                                        double mu_rate = 0.003;
                                        int thresold = Thresold_s[Tind];
                                        double regulation_parameter = Reg_Par[RP];
                                        //double mu_rate = muindex[Convert.ToInt32(MM)];


                                        //ss2 = "fit_N_" + N.ToString() + "pop_" + NPop.ToString() + method + "_cr_" + cr_rate.ToString() + "_mu_rate_" + mu_rate.ToString() + "Tre_" + thresold.ToString() + "_Reg_par" + regulation_parameter.ToString() + "_test_" + Ntest.ToString() + "OF_" + OF.ToString() + ".txt";
                                        #region file_existency
                                        if (!File.Exists(result))
                                        {
                                            for (int tt = 0; tt < Ntest; tt++)
                                            {
                                                NF[] GA = new NF[NPop];
                                                NF Best_indi = new NF();
                                                Best_indi.F = -1000000000;
                                                Fmax = -100000000; Fmin = 10000000000;
                                                ////////////////////////////////////
                                                for (int ii = 0; ii < 3; ii++)
                                                {
                                                    SP[ii] = 1.0 / 3;
                                                    Crossover_methods[ii] = 1.0 / 3;
                                                }
                                                //////////////////////////////
                                                for (int ii = 0; ii < NPop; ii++)
                                                {
                                                    GA[ii].x = NFC.node_constructing(N, Rand);
                                                    if (OF < 16 || OF == 22)
                                                        GA[ii].F = NFC.compute_fitness(GA[ii].x, OF, N);
                                                    else
                                                        GA[ii].F = NFC.compute_fitness(GA[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                    if (Fmax < GA[ii].F)
                                                        Fmax = GA[ii].F;
                                                }
                                                for (int it = 0; it < Iteration; it++)
                                                {
                                                    Fmax = -1000000000; Fmin = 100000000;
                                                    ////////////////////////////////////////////////////
                                                    for (int ii = 0; ii < NPop; ii++)
                                                    {

                                                        if (OF < 16 || OF == 22)
                                                            GA[ii].F = NFC.compute_fitness(GA[ii].x, OF, N);
                                                        else
                                                            GA[ii].F = NFC.compute_fitness(GA[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                        if (Fmax < GA[ii].F)
                                                            Fmax = GA[ii].F;
                                                        if (Fmin > GA[ii].F)
                                                            Fmin = GA[ii].F;
                                                    }
                                                    if (Fmin == Fmax)
                                                    {
                                                        for (int ii = 0; ii < NPop; ii++)
                                                        {
                                                            GA[ii].x = NFC.node_constructing(N, Rand);
                                                            if (OF < 16 || OF == 22)
                                                                GA[ii].F = NFC.compute_fitness(GA[ii].x, OF, N);
                                                            else
                                                                GA[ii].F = NFC.compute_fitness(GA[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                            if (Fmax < GA[ii].F)
                                                                Fmax = GA[ii].F;
                                                            if (Fmin > GA[ii].F)
                                                                Fmin = GA[ii].F;
                                                        }
                                                    }
                                                    if (Fmax > Best_indi.F)
                                                    {
                                                        Best_indi.F = Fmax;
                                                        Best_indi.x = (double[])GA[min_index].x.Clone();
                                                    }
                                                    int index_Selection = 0; int index_cross_over = 0;
                                                    /////////////Selecting the best method for selection operator
                                                    if (it < 1)
                                                    {
                                                        index_Selection = Rand.Next(SP.Length);
                                                        index_Selection = Rand.Next(Crossover_methods.Length);
                                                    }
                                                    else
                                                    {
                                                        double max_SP = -100; double max_Crov = -100;
                                                        for (int index = 0; index < SP.Length; index++)
                                                        {
                                                            if (max_SP < SP[index])
                                                            {
                                                                max_SP = SP[index];
                                                                index_Selection = index;
                                                            }
                                                            if (max_Crov < Crossover_methods[index])
                                                            {
                                                                max_Crov = Crossover_methods[index];
                                                                index_cross_over = index;
                                                            }
                                                        }
                                                    }
                                                    double[] OldCost = new double[NPop];
                                                    for (int ii = 0; ii < NPop; ii++)
                                                    {
                                                        OldCost[ii] = GA[ii].F;
                                                        if (cr_rate > Rand.NextDouble())
                                                        {
                                                            int index = NFC.Selection_Operators(index_Selection, GA, Fmin, NPop, Rand);
                                                            int index2 = NFC.Selection_Operators(index_Selection, GA, Fmin, NPop, Rand);
                                                            NF childs = NFC.Crossover_Operators(GA[index], GA[index2], N, OF, Rand, index_cross_over);
                                                            GA[ii].x = (double[])childs.x.Clone();
                                                        }
                                                        NFC.mutation(ref GA[ii], mu_rate, N, Rand);
                                                        if (OF < 16 || OF == 22)
                                                            GA[ii].F = NFC.compute_fitness(GA[ii].x, OF, N);
                                                        else
                                                            GA[ii].F = NFC.compute_fitness(GA[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                    }
                                                    /////Update_Operator////////////////////
                                                    double mean = NFC.update_operator_SRE(GA, OldCost, thresold);
                                                    if (mean > thresold)
                                                    {
                                                        SP[index_Selection] /= regulation_parameter;
                                                        Crossover_methods[index_cross_over] /= regulation_parameter;
                                                    }
                                                    else
                                                    {
                                                        SP[index_Selection] *= regulation_parameter;
                                                        Crossover_methods[index_cross_over] *= regulation_parameter;
                                                    }
                                                    ///////////////////////////////////////////
                                                }
                                                costs[tt] = Best_indi.F;
                                                Nodes_test[tt].x = (double[])Best_indi.x.Clone();
                                                Range.Text = mu_rate.ToString() + " " + cr_rate.ToString() + " " + OF.ToString() + " " + tt.ToString();
                                                Range.Refresh();
                                            }

                                            sw = new StreamWriter(result);
                                            for (int tt = 0; tt < Ntest; tt++)
                                            {
                                                sw.Write(costs[tt].ToString() + " ");
                                            }
                                            sw.Close();

                                        }
                                        #endregion
                                    }

                            }
                            #endregion
                            #region MASSWC
                            else if (method == "MASSWC")
                            {
                                 O = NPop * Convert.ToInt32(NI.Text);
                                //double[] N_LS_GL = { 0.1, 0.3, 0.5, 0.7, 0.9 };
                                //double[] Mu_rate_s = {0.00001,0.001,0.01,0.1,0.5 };
                                double[] N_LS_GL = { 0.7, 0.1, 0.5, 0.7, 0.7, 0.9 };
                                double[] Mu_rate_s = { 0.05, 0.0, 0.0, 0.0, 0.01, 0.1 };
                                //string[] Structs = { "Ring", "Star", "Btree", "Bipartie", "Cluster", "Grid", "Ladder", "Cross-Ladder", "Cellular" };
                                string[] Structs = { "Cellular" };
                                //for (int Portion_index = 0; Portion_index < 5; Portion_index++)
                                //    for (int M_index = 0; M_index < 5; M_index++)
                                //    {
                                for (int SS = 0; SS < 1; SS++)
                                {
                                    string Struct = Structs[SS];
                                    //double[] N_LS_GL = { 0.3,0.5,0.5,0.3,0.5,0.5,0.5,0.5,0.3,0.5,0.7,0.3,0.5,0,0.3,0.5,0.5,0.5,0.5,0.7,0.3,0.5 };
                                    //double[] Mu_rate_s = { 5, 0.5, 0.001, 0.001, 5, 5, 0.5, 0.5, 0.001, 0.01, 0.001, 0.01, 0.5, 0, 0.001, 0.1, 0.1, 0.5, 0.5, 5, 5, 0.1 };
                                    double LS_LG = N_LS_GL[OF - 22];
                                    double Mu_rate = Mu_rate_s[OF - 22];
                                    Iteration = Convert.ToInt32((O * (1 - LS_LG)) / NPop) + 1;
                                    int I_str = Convert.ToInt32(O * LS_LG / Iteration);
                                    int best_index = 0;
                                    bool local_search = false;
                                    //if (Mu_rate == 0.00001)
                                    //    ss2 = "fit_N_" + N.ToString() + "_Mr_" + 5.ToString() + "LSLG" + LS_LG.ToString() + "_Pop_" + NPop.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "Struct_" + Struct + "Lap" + "_OF" + OF.ToString() + ".txt";
                                    //else
                                    //    ss2 = "fit_N_" + N.ToString() + "_Mr_" + Mu_rate.ToString() + "LSLG" + LS_LG.ToString() + "_Pop_" + NPop.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "Struct_" + Struct + "Lap" + "_OF" + OF.ToString() + ".txt";
                                    if (!File.Exists(result))
                                    {
                                        for (int tt = 0; tt < Ntest; tt++)
                                        {
                                            LS[] Pop = new LS[NPop];
                                            NF Best_nodes = new NF();
                                            NF[] Breed = new NF[NPop];
                                            Best_nodes.F = -1E250;
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                Pop[ii].x = NFC.node_constructing(N, Rand);

                                                if (OF < 16 || OF == 22)
                                                    Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N);
                                                else
                                                    Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                Pop[ii].local_search = true;

                                            }
                                            for (int it = 0; it < Iteration; it++)
                                            {
                                                //perform Steady State GA
                                                #region generate Breeds
                                                double max = -1E250; int max_index = 0;
                                                local_search = false;
                                                for (int ii = 0; ii < NPop; ii++)
                                                {
                                                    int parent1 = ii;
                                                    int[] S = NFC.Structures(ii, NPop, Struct);
                                                    int index2 = Rand.Next(S.Length);
                                                    int parent2 = S[index2];
                                                    Breed[ii].x = NFC.BLX_Crossover(Pop[parent1].x, Pop[parent2].x);
                                                    if (Rand.NextDouble() < Mu_rate)
                                                        NFC.Mutation_BGA2(ref Breed[ii], Pop, N, Rand);
                                                    if (OF < 16 || OF == 22)
                                                        Breed[ii].F = NFC.compute_fitness(Breed[ii].x, OF, N);
                                                    else
                                                        Breed[ii].F = NFC.compute_fitness(Breed[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                    Breed[ii].x = NFC.Node_treatment(Breed[ii].x, N, Rand);
                                                    if (Pop[ii].F < Breed[ii].F)
                                                    {
                                                        Pop[ii].x = (double[])Breed[ii].x.Clone();
                                                        Pop[ii].F = Breed[ii].F;
                                                        Pop[ii].visited = 0;
                                                        Pop[ii].local_search = true;
                                                    }
                                                    if (Pop[ii].F > max && Pop[ii].local_search == true)
                                                    {
                                                        max = Pop[ii].F;
                                                        max_index = ii;
                                                        local_search = true;
                                                    }

                                                }
                                                #endregion
                                                // Replacement Strategy

                                                //
                                                //Pick the best individual of Population
                                                if (local_search == false)
                                                {
                                                    for (int ii = 1; ii < NPop; ii++)
                                                    {
                                                        Pop[ii].x = NFC.node_constructing(N, Rand);
                                                        if (OF < 16 || OF == 22)
                                                            Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N);
                                                        else
                                                            Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                        Pop[ii].local_search = true;
                                                        Pop[ii].visited = 0;
                                                        if (Pop[ii].F > max && Pop[ii].local_search == true)
                                                        {
                                                            max = Pop[ii].F;
                                                            max_index = ii;
                                                        }
                                                    }
                                                }
                                                double diversity = NFC.diversity(Pop, NPop, N);
                                                if (diversity < 0.01)
                                                {
                                                    for (int ii = 1; ii < NPop; ii++)
                                                    {
                                                        Pop[ii].x = NFC.node_constructing(N, Rand);
                                                        if (OF < 16 || OF == 22)
                                                            Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N);
                                                        else
                                                            Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                        Pop[ii].local_search = true;
                                                        Pop[ii].visited = 0;
                                                        if (Pop[ii].F > max && Pop[ii].local_search == true)
                                                        {
                                                            max = Pop[ii].F;
                                                            max_index = ii;
                                                        }
                                                    }
                                                }
                                                LS Best_Indi_BLS = new LS();
                                                Best_Indi_BLS.F = Pop[max_index].F;
                                                Best_Indi_BLS.x = (double[])Pop[max_index].x.Clone();
                                                Best_Indi_BLS.visited = Pop[max_index].visited;
                                                Best_Indi_BLS.rho = Pop[max_index].rho;
                                                if (Pop[max_index].visited == 1)
                                                    Best_Indi_BLS.bias = (double[])Pop[max_index].bias.Clone();
                                                best_index = max_index;
                                                // Initialize Bias and Rho
                                                if (Best_Indi_BLS.visited == 0)
                                                {
                                                    Best_Indi_BLS.bias = new double[N];
                                                    Best_Indi_BLS.rho = 0.01;
                                                }
                                                // Local Search Operator
                                                LS Best_Indi_ALS = new LS();
                                                Best_Indi_ALS.x = (double[])Best_Indi_BLS.x.Clone();
                                                Best_Indi_ALS.F = Best_Indi_BLS.F;
                                                Best_Indi_ALS.rho = Best_Indi_BLS.rho;
                                                Best_Indi_ALS.visited = Best_Indi_BLS.visited;
                                                Best_Indi_ALS.bias = (double[])Best_Indi_BLS.bias.Clone();
                                                NFC.SolisWets(ref Best_Indi_ALS.x, ref Best_Indi_ALS.F, N, ref Best_Indi_ALS.bias, ref Best_Indi_ALS.rho, I_str, Rand, OF, Perm, Orthogonal_Matrix, M, o);

                                                // Remove Unchanged particle from population
                                                if (Best_Indi_BLS.F - Best_Indi_ALS.F == 0)
                                                {
                                                    Pop[best_index].local_search = false;
                                                    //  Pop[best_index].visited = 1;
                                                }
                                                else
                                                {
                                                    Pop[best_index].x = (double[])Best_Indi_ALS.x.Clone();
                                                    Pop[best_index].visited = 1;
                                                    Pop[best_index].F = Best_Indi_ALS.F;
                                                    Pop[best_index].rho = Best_Indi_ALS.rho;
                                                    Pop[best_index].bias = (double[])Best_Indi_ALS.bias.Clone();
                                                    Pop[best_index].local_search = true;
                                                }

                                                Pop = NFC.mergesort(Pop, NPop);
                                                if (Best_nodes.F < Pop[0].F)
                                                {
                                                    Best_nodes.F = Pop[0].F;
                                                    Best_nodes.x = (double[])Pop[0].x.Clone();
                                                }

                                                /////////////////////

                                            }
                                            costs[tt] = Best_nodes.F;

                                        }


                                        sw = new StreamWriter(result);
                                        for (int tt = 0; tt < Ntest; tt++)
                                        {
                                            sw.Write(costs[tt].ToString() + " ");
                                        }
                                        sw.Close();
                                    }
                                }
                                //}

                            }

                            #endregion
                            #region 3LMA
                            else if (method == "3TMA")
                            {

                                //function_evaluation2013 class_new = new function_evaluation2013(ub, lb, dimension);
                                //double O = NPop * Convert.ToInt32(NI.Text);
                                //double[] N_LS_GL = { 0.3, 0.1, 0.3, 0.9, 0.5, 0.5, 0.7, 0.3, 0.9, 0.1, 0.3, 0.9, 0.7, 0, 0.9, 0.3, 0.7, 0.5, 0.9, 0.3, 0.9, 0.3 };
                                //double[] Mu_rate_s = { 0.1, 0.8, 0.5, 0.8, 0.8, 1, 1, 1, 0.8, 1, 0.8, 0.8, 0.8, 0, 1, 1, 0.8, 0.8, 1, 0.5, 1, 0.01 };
                                ////double[] N_LS_GL = { 0.5, 0.7 };
                                ////double[] Mu_rate_s = { 0.01, 0.1, 0.5 };
                                //double[] Jrs = { 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };
                                ////string[] Structs = { "Ring", "Star", "Btree", "Bipartie", "Cluster", "Grid", "Ladder", "Cross-Ladder", "Cellular" };
                                //string[] Structs = { "Star" };
                                ////for (int Portion_index = 0; Portion_index < 5; Portion_index++)
                                ////    for (int M_index = 0; M_index < 1; M_index++)
                                ////    {

                                //bool OBL = true;
                                ////for (int OB = 0; OB < 10; OB++)
                                ////{
                                //for (int SS = 0; SS < 1; SS++)
                                //{
                                //    string Struct = Structs[SS];
                                //    //double[] N_LS_GL = { 0.3, 0.5, 0.9, 0.5, 0.5, 0.3, 0.3, 0.5, 0.5, 0.3, 0.7, 0.5, 0.5, 0, 0.5, 0.7, 0.9, 0.9, 0.9, 0.7, 0.5, 0.7 };
                                //    //double[] Mu_rate_s = { 0.1, 0.00001, 5, 0.1, 0.1, 0.5, 0.01, 0.1, 5, 0.5, 0.001, 5, 0.1, 0, 0.1, 5, 0.01, 0.5, 0.01, 0.1, 0.5, 0.5 };
                                //    double LS_LG = N_LS_GL[OF - 1];
                                //    double Mu_rate = Mu_rate_s[OF - 1];
                                //    Iteration = Convert.ToInt32((O * (1 - LS_LG)) / NPop) + 1;
                                //    //int I_str = Convert.ToInt32(O * LS_LG) / Iteration;
                                //    int I_str = 100;
                                //    int best_index = 0;
                                //    bool local_search = false;
                                //    double Jr = 0.3;
                                //    ss2 = "fit_N_" + N.ToString() + "_Mr_" + Mu_rate.ToString() + "LSLG" + LS_LG.ToString() + "_J" + Jr.ToString() + "_Pop_" + NPop.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "Struct_" + Struct + "_OF" + OF.ToString() + ".txt";
                                //    if (!File.Exists(ss2))
                                //    {
                                //        for (int tt = 0; tt < Ntest; tt++)
                                //        {

                                //            //counter_JR = 1;
                                //            LS[] Pop = new LS[NPop];
                                //            LS[] OPop = new LS[NPop];
                                //            NF Best_nodes = new NF();
                                //            NF[] Breed = new NF[NPop];
                                //            Best_nodes.F = 1E250;
                                //            for (int ii = 0; ii < NPop; ii++)
                                //            {
                                //                Pop[ii].x = NFC.node_constructing2(N, Rand, ub, lb);
                                //                if (OBL)
                                //                    OPop[ii].x = NFC.Opposition_Learning(Pop[ii].x);
                                //                Pop[ii].F = class_new.Choosing_Function(OF, Pop[ii].x, dimension, xopt, Perm, r25, r50, r100, S_Array, W_Array, overlap);
                                //                if (OBL)
                                //                {
                                //                    OPop[ii].F = class_new.Choosing_Function(OF, OPop[ii].x, dimension, xopt, Perm, r25, r50, r100, S_Array, W_Array, overlap);
                                //                    if (OPop[ii].F < Pop[ii].F)
                                //                    {
                                //                        Pop[ii].F = OPop[ii].F;
                                //                        Pop[ii].x = (double[])OPop[ii].x.Clone();
                                //                    }
                                //                }
                                //                Pop[ii].local_search = true;

                                //            }
                                //            int OBL_iter = 0;
                                //            for (int it = 0; it < Iteration; it++)
                                //            {
                                //                //perform Steady State GA
                                //                #region generate Breeds
                                //                double min = 1E250; int min_index = 0;
                                //                local_search = false;
                                //                int imax = -1; int imin = -1; double max_b = -1E250; double min_p = 1E250;
                                //                for (int ii = 0; ii < NPop; ii++)
                                //                {
                                //                    int parent1 = Rand.Next(NPop);
                                //                    int[] S = NFC.Structures(ii, NPop, Struct);
                                //                    int index2 = Rand.Next(S.Length);
                                //                    //int parent2 = S[index2];
                                //                    int parent2 = NFC.negative_assortative_mating_strategy(Pop[parent1].x, Pop, Rand, S);
                                //                    Breed[ii].x = NFC.BLX_Crossover(Pop[parent1].x, Pop[parent2].x);
                                //                    if (Rand.NextDouble() < Mu_rate)
                                //                        NFC.Mutation_BGA2(ref Breed[ii], Pop, N, Rand);
                                //                    Breed[ii].F = class_new.Choosing_Function(OF, Pop[ii].x, dimension, xopt, Perm, r25, r50, r100, S_Array, W_Array, overlap);
                                //                    Breed[ii].x = NFC.Node_treatment_new(Breed[ii].x, N, Rand, 100, -100);
                                //                    ////////////////////////////////////////////
                                //                    if (Breed[ii].F < min_p)
                                //                    {
                                //                        imin = ii;
                                //                        min_p = Breed[ii].F;
                                //                    }
                                //                    if (Pop[ii].F > max_b)
                                //                    {
                                //                        imax = ii;
                                //                        max_b = Pop[ii].F;
                                //                    }
                                //                }
                                //                if (max_b > min_p)
                                //                {
                                //                    Pop[imax].F = Breed[imin].F;
                                //                    Pop[imax].x = (double[])Breed[imin].x.Clone();
                                //                }
                                //                /////////////////////////////**** The second layer of the algorithm
                                //                for (int ii = 0; ii < NPop; ii++)
                                //                {
                                //                    if (OBL)
                                //                    {
                                //                        if (Jr > Rand.NextDouble())
                                //                        {

                                //                            OPop[ii].x = NFC.generating_jumping(Pop, Pop[ii].x, NPop, N, ii);
                                //                            OPop[ii].F = class_new.Choosing_Function(OF, OPop[ii].x, dimension, xopt, Perm, r25, r50, r100, S_Array, W_Array, overlap);
                                //                            if (Pop[ii].F > OPop[ii].F)
                                //                            {
                                //                                Pop[ii].x = (double[])OPop[ii].x.Clone();
                                //                                Pop[ii].F = OPop[ii].F;
                                //                            }
                                //                            OBL_iter++;
                                //                            if (OBL_iter == NPop)
                                //                            {
                                //                                OBL_iter = 0;
                                //                                it++;
                                //                            }
                                //                        }
                                //                    }
                                //                    /////////////////////////////////////////////////////////////////
                                //                    if (Pop[ii].F < min && Pop[ii].local_search == true)
                                //                    {
                                //                        min = Pop[ii].F;
                                //                        min_index = ii;
                                //                        local_search = true;
                                //                    }
                                //                }
                                //                #endregion
                                //                // Replacement Strategy

                                //                //
                                //                //Pick the best individual of Population
                                //                if (local_search == false)
                                //                {

                                //                    for (int ii = 0; ii < NPop; ii++)
                                //                    {
                                //                        Pop[ii].x = NFC.node_constructing2(N, Rand, 100, -100);
                                //                        Pop[ii].F = class_new.Choosing_Function(OF, Pop[ii].x, dimension, xopt, Perm, r25, r50, r100, S_Array, W_Array, overlap);
                                //                        Pop[ii].x = NFC.Node_treatment(Pop[ii].x, N, Rand);
                                //                        Pop[ii].local_search = true;
                                //                        Pop[ii].visited = 0;
                                //                        if (Pop[ii].F < min && Pop[ii].local_search == true)
                                //                        {
                                //                            min = Pop[ii].F;
                                //                            min_index = ii;
                                //                        }
                                //                    }
                                //                    int randx = Rand.Next(0, NPop);
                                //                    Pop[randx].x = (double[])Best_nodes.x.Clone();
                                //                    Pop[randx].F = Best_nodes.F;
                                //                }
                                //                LS Best_Indi_BLS = new LS();
                                //                Best_Indi_BLS.F = Pop[min_index].F;
                                //                Best_Indi_BLS.x = (double[])Pop[min_index].x.Clone();
                                //                Best_Indi_BLS.visited = Pop[min_index].visited;
                                //                Best_Indi_BLS.rho = Pop[min_index].rho;
                                //                if (Pop[min_index].visited == 1)
                                //                    Best_Indi_BLS.bias = (double[])Pop[min_index].bias.Clone();
                                //                best_index = min_index;
                                //                // Initialize Bias and Rho
                                //                if (Best_Indi_BLS.visited == 0)
                                //                {
                                //                    Best_Indi_BLS.bias = new double[N];
                                //                    Best_Indi_BLS.rho = 0.01;
                                //                }
                                //                // Local Search Operator
                                //                LS Best_Indi_ALS = new LS();
                                //                Best_Indi_ALS.x = (double[])Best_Indi_BLS.x.Clone();
                                //                Best_Indi_ALS.F = Best_Indi_BLS.F;
                                //                Best_Indi_ALS.rho = Best_Indi_BLS.rho;
                                //                Best_Indi_ALS.visited = Best_Indi_BLS.visited;
                                //                Best_Indi_ALS.bias = (double[])Best_Indi_BLS.bias.Clone();
                                //                //if(switches==0)
                                //                NFC.SolisWets(ref Best_Indi_ALS.x, ref Best_Indi_ALS.F, N, ref Best_Indi_ALS.bias, ref Best_Indi_ALS.rho, I_str, Rand, OF, xopt, dimension, Perm, r25, r50, r100, S_Array, W_Array, overlap, ub, lb);
                                //                //else
                                //                //    NFC.Subgrouping_SolisWets(ref Best_Indi_ALS.x, ref Best_Indi_ALS.F, N, ref Best_Indi_ALS.bias, ref Best_Indi_ALS.rho, I_str, Rand, OF, Perm, Orthogonal_Matrix, M, o);
                                //                // Remove Unchanged particle from population
                                //                if (Best_Indi_BLS.F - Best_Indi_ALS.F == 0)
                                //                {
                                //                    Pop[best_index].local_search = false;
                                //                    //  Pop[best_index].visited = 1;
                                //                }
                                //                else
                                //                {
                                //                    Pop[best_index].x = (double[])Best_Indi_ALS.x.Clone();
                                //                    Pop[best_index].visited = 1;
                                //                    Pop[best_index].F = Best_Indi_ALS.F;
                                //                    Pop[best_index].rho = Best_Indi_ALS.rho;
                                //                    Pop[best_index].bias = (double[])Best_Indi_ALS.bias.Clone();
                                //                    Pop[best_index].local_search = true;
                                //                }

                                //                if (Best_nodes.F > Pop[best_index].F)
                                //                {
                                //                    Best_nodes.F = Pop[best_index].F;
                                //                    Best_nodes.x = (double[])Pop[best_index].x.Clone();
                                //                }

                                //                /////////////////////
                                //                //counter_JR=counter_JR * 0.99;
                                //                //Jr = 1 - (counter_JR);
                                //            }
                                //            costs[tt] = Best_nodes.F;

                                //        }


                                //        sw = new StreamWriter(ss2);
                                //        for (int tt = 0; tt < Ntest; tt++)
                                //        {
                                //            sw.Write(costs[tt].ToString() + " ");
                                //        }
                                //        sw.Close();
                                //}
                                //}
                                //}
                            }

                            #endregion
                            #region 3LMA_OLD
                            else if (method == "3LMA_OLD")
                            {
                                O = NPop * Convert.ToInt32(NI.Text);
                                double[] N_LS_GL = { 0.3, 0.1, 0.3, 0.9, 0.5, 0.5, 0.7, 0.3, 0.9, 0.1, 0.3, 0.9, 0.7, 0, 0.9, 0.3, 0.7, 0.5, 0.9, 0.3, 0.9, 0.3 };
                                double[] Mu_rate_s = { 0.1, 0.8, 0.5, 0.8, 0.8, 1, 1, 1, 0.8, 1, 0.8, 0.8, 0.8, 0, 1, 1, 0.8, 0.8, 1, 0.5, 1, 0.01 };
                                //double[] N_LS_GL = { 0.5, 0.7 };
                                //double[] Mu_rate_s = { 0.01, 0.1, 0.5 };
                                double[] Jrs = { 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };
                                //string[] Structs = { "Ring", "Star", "Btree", "Bipartie", "Cluster", "Grid", "Ladder", "Cross-Ladder", "Cellular" };
                                string[] Structs = { "Star" };
                                //for (int Portion_index = 0; Portion_index < 5; Portion_index++)
                                //    for (int M_index = 0; M_index < 1; M_index++)
                                //    {

                                bool OBL = true;
                                //for (int OB = 0; OB < 10; OB++)
                                //{
                                for (int SS = 0; SS < 1; SS++)
                                {
                                    string Struct = Structs[SS];
                                    //double[] N_LS_GL = { 0.3, 0.5, 0.9, 0.5, 0.5, 0.3, 0.3, 0.5, 0.5, 0.3, 0.7, 0.5, 0.5, 0, 0.5, 0.7, 0.9, 0.9, 0.9, 0.7, 0.5, 0.7 };
                                    //double[] Mu_rate_s = { 0.1, 0.00001, 5, 0.1, 0.1, 0.5, 0.01, 0.1, 5, 0.5, 0.001, 5, 0.1, 0, 0.1, 5, 0.01, 0.5, 0.01, 0.1, 0.5, 0.5 };
                                    double LS_LG = N_LS_GL[OF - 1];
                                    double Mu_rate = Mu_rate_s[OF - 1];
                                    Iteration = Convert.ToInt32((O * (1 - LS_LG)) / NPop) + 1;
                                    int I_str = Convert.ToInt32(O * LS_LG) / Iteration;
                                    int best_index = 0;
                                    bool local_search = false;
                                    double Jr = 0.3;
                                    //ss2 = "fit_N_" + N.ToString() + "_Mr_" + Mu_rate.ToString() + "LSLG" + LS_LG.ToString() + "_J" + Jr.ToString() + "_Pop_" + NPop.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "Struct_" + Struct + "_OF" + OF.ToString() + ".txt";
                                    if (!File.Exists(result))
                                    {
                                        for (int tt = 0; tt < Ntest; tt++)
                                        {

                                            //counter_JR = 1;
                                            LS[] Pop = new LS[NPop];
                                            LS[] OPop = new LS[NPop];
                                            NF Best_nodes = new NF();
                                            NF[] Breed = new NF[NPop];
                                            Best_nodes.F = -1E250;
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                Pop[ii].x = NFC.node_constructing(N, Rand);
                                                if (OBL)
                                                    OPop[ii].x = NFC.Opposition_Learning(Pop[ii].x);
                                                if (OF < 16 || OF == 22)
                                                    Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N);
                                                else
                                                    Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                if (OBL)
                                                {
                                                    if (OF < 16 || OF == 22)
                                                        OPop[ii].F = NFC.compute_fitness(OPop[ii].x, OF, N);
                                                    else
                                                        OPop[ii].F = NFC.compute_fitness(OPop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);

                                                    if (OPop[ii].F > Pop[ii].F)
                                                    {
                                                        Pop[ii].F = OPop[ii].F;
                                                        Pop[ii].x = (double[])OPop[ii].x.Clone();
                                                    }
                                                }
                                                Pop[ii].local_search = true;

                                            }
                                            int OBL_iter = 0;
                                            for (int it = 0; it < Iteration; it++)
                                            {
                                                //perform Steady State GA
                                                #region generate Breeds
                                                double max = -1E250; int max_index = 0;
                                                local_search = false;
                                                int imax = -1; int imin = -1; double max_b = -1E250; double min_p = 1E250;
                                                for (int ii = 0; ii < NPop; ii++)
                                                {
                                                    int parent1 = Rand.Next(NPop);
                                                    int[] S = NFC.Structures(ii, NPop, Struct);
                                                    int index2 = Rand.Next(S.Length);
                                                    //int parent2 = S[index2];
                                                    int parent2 = NFC.negative_assortative_mating_strategy(Pop[parent1].x, Pop, Rand, S);
                                                    Breed[ii].x = NFC.BLX_Crossover(Pop[parent1].x, Pop[parent2].x);
                                                    if (Rand.NextDouble() < Mu_rate)
                                                        NFC.Mutation_BGA2(ref Breed[ii], Pop, N, Rand);
                                                    if (OF < 16 || OF == 22)
                                                        Breed[ii].F = NFC.compute_fitness(Breed[ii].x, OF, N);
                                                    else
                                                        Breed[ii].F = NFC.compute_fitness(Breed[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                    Breed[ii].x = NFC.Node_treatment(Breed[ii].x, N, Rand);
                                                    ////////////////////////////////////////////
                                                    if (Breed[ii].F > max_b)
                                                    {
                                                        imax = ii;
                                                        max_b = Breed[ii].F;
                                                    }
                                                    if (Pop[ii].F < min_p)
                                                    {
                                                        imin = ii;
                                                        min_p = Pop[ii].F;
                                                    }
                                                }
                                                if (max_b > min_p)
                                                {
                                                    Pop[imin].F = Breed[imax].F;
                                                    Pop[imin].x = (double[])Breed[imax].x.Clone();
                                                }
                                                /////////////////////////////**** The second layer of the algorithm
                                                for (int ii = 0; ii < NPop; ii++)
                                                {
                                                    if (OBL)
                                                    {
                                                        if (Jr > Rand.NextDouble())
                                                        {

                                                            OPop[ii].x = NFC.generating_jumping(Pop, Pop[ii].x, NPop, N, ii);
                                                            if (OF < 16 || OF == 22)
                                                                OPop[ii].F = NFC.compute_fitness(OPop[ii].x, OF, N);
                                                            else
                                                                OPop[ii].F = NFC.compute_fitness(OPop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                            if (Pop[ii].F < OPop[ii].F)
                                                            {
                                                                Pop[ii].x = (double[])OPop[ii].x.Clone();
                                                                Pop[ii].F = OPop[ii].F;
                                                            }
                                                            OBL_iter++;
                                                            if (OBL_iter == NPop)
                                                            {
                                                                OBL_iter = 0;
                                                                it++;
                                                            }
                                                        }
                                                    }
                                                    /////////////////////////////////////////////////////////////////
                                                    if (Pop[ii].F > max && Pop[ii].local_search == true)
                                                    {
                                                        max = Pop[ii].F;
                                                        max_index = ii;
                                                        local_search = true;
                                                    }
                                                }
                                                #endregion
                                                // Replacement Strategy

                                                //
                                                //Pick the best individual of Population
                                                if (local_search == false)
                                                {

                                                    for (int ii = 0; ii < NPop; ii++)
                                                    {
                                                        Pop[ii].x = NFC.node_constructing(N, Rand);
                                                        if (OF < 16 || OF == 22)
                                                            Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N);
                                                        else
                                                            Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                        Pop[ii].x = NFC.Node_treatment(Pop[ii].x, N, Rand);
                                                        Pop[ii].local_search = true;
                                                        Pop[ii].visited = 0;
                                                        if (Pop[ii].F > max && Pop[ii].local_search == true)
                                                        {
                                                            max = Pop[ii].F;
                                                            max_index = ii;
                                                        }
                                                    }
                                                    int randx = Rand.Next(0, NPop);
                                                    Pop[randx].x = (double[])Best_nodes.x.Clone();
                                                    Pop[randx].F = Best_nodes.F;
                                                }
                                                LS Best_Indi_BLS = new LS();
                                                Best_Indi_BLS.F = Pop[max_index].F;
                                                Best_Indi_BLS.x = (double[])Pop[max_index].x.Clone();
                                                Best_Indi_BLS.visited = Pop[max_index].visited;
                                                Best_Indi_BLS.rho = Pop[max_index].rho;
                                                if (Pop[max_index].visited == 1)
                                                    Best_Indi_BLS.bias = (double[])Pop[max_index].bias.Clone();
                                                best_index = max_index;
                                                // Initialize Bias and Rho
                                                if (Best_Indi_BLS.visited == 0)
                                                {
                                                    Best_Indi_BLS.bias = new double[N];
                                                    Best_Indi_BLS.rho = 0.01;
                                                }
                                                // Local Search Operator
                                                LS Best_Indi_ALS = new LS();
                                                Best_Indi_ALS.x = (double[])Best_Indi_BLS.x.Clone();
                                                Best_Indi_ALS.F = Best_Indi_BLS.F;
                                                Best_Indi_ALS.rho = Best_Indi_BLS.rho;
                                                Best_Indi_ALS.visited = Best_Indi_BLS.visited;
                                                Best_Indi_ALS.bias = (double[])Best_Indi_BLS.bias.Clone();
                                                //if(switches==0)
                                                NFC.SolisWets(ref Best_Indi_ALS.x, ref Best_Indi_ALS.F, N, ref Best_Indi_ALS.bias, ref Best_Indi_ALS.rho, I_str, Rand, OF, Perm, Orthogonal_Matrix, M, o);
                                                //else
                                                //    NFC.Subgrouping_SolisWets(ref Best_Indi_ALS.x, ref Best_Indi_ALS.F, N, ref Best_Indi_ALS.bias, ref Best_Indi_ALS.rho, I_str, Rand, OF, Perm, Orthogonal_Matrix, M, o);
                                                // Remove Unchanged particle from population
                                                if (Best_Indi_BLS.F - Best_Indi_ALS.F == 0)
                                                {
                                                    Pop[best_index].local_search = false;
                                                    //  Pop[best_index].visited = 1;
                                                }
                                                else
                                                {
                                                    Pop[best_index].x = (double[])Best_Indi_ALS.x.Clone();
                                                    Pop[best_index].visited = 1;
                                                    Pop[best_index].F = Best_Indi_ALS.F;
                                                    Pop[best_index].rho = Best_Indi_ALS.rho;
                                                    Pop[best_index].bias = (double[])Best_Indi_ALS.bias.Clone();
                                                    Pop[best_index].local_search = true;
                                                }

                                                if (Best_nodes.F < Pop[best_index].F)
                                                {
                                                    Best_nodes.F = Pop[best_index].F;
                                                    Best_nodes.x = (double[])Pop[best_index].x.Clone();
                                                }

                                                /////////////////////
                                                //counter_JR=counter_JR * 0.99;
                                                //Jr = 1 - (counter_JR);
                                            }
                                            costs[tt] = Best_nodes.F;

                                        }


                                        sw = new StreamWriter(result);
                                        for (int tt = 0; tt < Ntest; tt++)
                                        {
                                            sw.Write(costs[tt].ToString() + " ");
                                        }
                                        sw.Close();
                                    }
                                    //}
                                }
                            }
                            #endregion
                            #region (1+1)LMA
                            else if (method == "(1+1)LMA")
                            {
                                O = NPop * Convert.ToInt32(NI.Text);
                                double[] N_LS_GL = { 0.3, 0.1, 0.3, 0.9, 0.5, 0.5, 0.7, 0.3, 0.9, 0.1, 0.3, 0.9, 0.7, 0, 0.9, 0.3, 0.7, 0.5, 0.9, 0.3, 0.9, 0.3 };
                                double[] Mu_rate_s = { 0.1, 0.8, 0.5, 0.8, 0.8, 1, 1, 1, 0.8, 1, 0.8, 0.8, 0.8, 0, 1, 1, 0.8, 0.8, 1, 0.5, 1, 0.01 };
                                //double[] N_LS_GL = { 0.5, 0.7 };
                                //double[] Mu_rate_s = { 0.01, 0.1, 0.5 };
                                // double[] Jrs = { 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };
                                //string[] Structs = { "Ring", "Star", "Btree", "Bipartie", "Cluster", "Grid", "Ladder", "Cross-Ladder", "Cellular" };
                                string[] Structs = { "Star" };
                                //for (int Portion_index = 0; Portion_index < 5; Portion_index++)
                                //    for (int M_index = 0; M_index < 1; M_index++)
                                //    {

                                bool OBL = false;
                                //for (int OB = 0; OB < 10; OB++)
                                //{
                                for (int SS = 0; SS < 1; SS++)
                                {
                                    string Struct = Structs[SS];
                                    //double[] N_LS_GL = { 0.3, 0.5, 0.9, 0.5, 0.5, 0.3, 0.3, 0.5, 0.5, 0.3, 0.7, 0.5, 0.5, 0, 0.5, 0.7, 0.9, 0.9, 0.9, 0.7, 0.5, 0.7 };
                                    //double[] Mu_rate_s = { 0.1, 0.00001, 5, 0.1, 0.1, 0.5, 0.01, 0.1, 5, 0.5, 0.001, 5, 0.1, 0, 0.1, 5, 0.01, 0.5, 0.01, 0.1, 0.5, 0.5 };
                                    double LS_LG = N_LS_GL[OF - 1];
                                    double Mu_rate = Mu_rate_s[OF - 1];
                                    Iteration = Convert.ToInt32((O * (1 - LS_LG)) / NPop) + 1;
                                    int I_str = Convert.ToInt32(O * LS_LG) / Iteration;
                                    int best_index = 0;
                                    bool local_search = false;
                                    double Jr = 0.3;
                                    double[,] Costs = new double[Ntest, Iteration / 2 + 1];
                                    //ss2 = "fit_N_" + N.ToString() + "_Mr_" + Mu_rate.ToString() + "LSLG" + LS_LG.ToString() + "_J" + Jr.ToString() + "_Pop_" + NPop.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "Struct_" + Struct + "_OF" + OF.ToString() + ".txt";
                                    if (!File.Exists(result))
                                    {
                                        for (int tt = 0; tt < Ntest; tt++)
                                        {

                                            //counter_JR = 1;
                                            LS[] Pop = new LS[NPop];
                                            LS[] OPop = new LS[NPop];
                                            NF Best_nodes = new NF();
                                            NF[] Breed = new NF[NPop];
                                            Best_nodes.F = -1E250;
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                Pop[ii].x = NFC.node_constructing(N, Rand);
                                                if (OBL)
                                                    OPop[ii].x = NFC.Opposition_Learning(Pop[ii].x);
                                                if (OF < 16 || OF == 22)
                                                    Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N);
                                                else
                                                    Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                if (OBL)
                                                {
                                                    if (OF < 16 || OF == 22)
                                                        OPop[ii].F = NFC.compute_fitness(OPop[ii].x, OF, N);
                                                    else
                                                        OPop[ii].F = NFC.compute_fitness(OPop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);

                                                    if (OPop[ii].F > Pop[ii].F)
                                                    {
                                                        Pop[ii].F = OPop[ii].F;
                                                        Pop[ii].x = (double[])OPop[ii].x.Clone();
                                                    }
                                                }
                                                Pop[ii].local_search = true;

                                            }

                                            int kk = 0;
                                            for (int it = 0; it < Iteration; it++)
                                            {
                                                //perform Steady State GA
                                                #region generate Breeds
                                                double max = -1E250; int max_index = 0;
                                                local_search = false;
                                                int imax = -1; int imin = -1; double max_b = -1E250; double min_p = 1E250;
                                                for (int ii = 0; ii < NPop; ii++)
                                                {
                                                    int parent1 = Rand.Next(NPop);
                                                    int[] S = NFC.Structures(ii, NPop, Struct);
                                                    int index2 = Rand.Next(S.Length);
                                                    //int parent2 = S[index2];
                                                    int parent2 = NFC.negative_assortative_mating_strategy(Pop[parent1].x, Pop, Rand, S);
                                                    Breed[ii].x = NFC.BLX_Crossover(Pop[parent1].x, Pop[parent2].x);
                                                    if (Rand.NextDouble() < Mu_rate)
                                                        NFC.Mutation_BGA2(ref Breed[ii], Pop, N, Rand);
                                                    if (OF < 16 || OF == 22)
                                                        Breed[ii].F = NFC.compute_fitness(Breed[ii].x, OF, N);
                                                    else
                                                        Breed[ii].F = NFC.compute_fitness(Breed[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                    Breed[ii].x = NFC.Node_treatment(Breed[ii].x, N, Rand);
                                                    ////////////////////////////////////////////
                                                    if (Breed[ii].F > max_b)
                                                    {
                                                        imax = ii;
                                                        max_b = Breed[ii].F;
                                                    }
                                                    if (Pop[ii].F < min_p)
                                                    {
                                                        imin = ii;
                                                        min_p = Pop[ii].F;
                                                    }
                                                }
                                                if (max_b > min_p)
                                                {
                                                    Pop[imin].F = Breed[imax].F;
                                                    Pop[imin].x = (double[])Breed[imax].x.Clone();
                                                }
                                                for (int ii = 0; ii < NPop; ii++)
                                                {
                                                    if (OBL)
                                                    {
                                                        if (Jr > Rand.NextDouble())
                                                        {

                                                            OPop[ii].x = NFC.generating_jumping(Pop, Pop[ii].x, NPop, N, ii);
                                                            if (OF < 16 || OF == 22)
                                                                OPop[ii].F = NFC.compute_fitness(OPop[ii].x, OF, N);
                                                            else
                                                                OPop[ii].F = NFC.compute_fitness(OPop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                            if (Pop[ii].F < OPop[ii].F)
                                                            {
                                                                Pop[ii].x = (double[])OPop[ii].x.Clone();
                                                                Pop[ii].F = OPop[ii].F;
                                                            }
                                                        }
                                                    }
                                                    if (Pop[ii].F > max && Pop[ii].local_search == true)
                                                    {
                                                        max = Pop[ii].F;
                                                        max_index = ii;
                                                        local_search = true;
                                                    }
                                                }
                                                #endregion
                                                // Replacement Strategy

                                                //
                                                //Pick the best individual of Population
                                                if (local_search == false)
                                                {

                                                    for (int ii = 0; ii < NPop; ii++)
                                                    {
                                                        Pop[ii].x = NFC.node_constructing(N, Rand);
                                                        if (OF < 16 || OF == 22)
                                                            Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N);
                                                        else
                                                            Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                        Pop[ii].x = NFC.Node_treatment(Pop[ii].x, N, Rand);
                                                        Pop[ii].local_search = true;
                                                        Pop[ii].visited = 0;
                                                        if (Pop[ii].F > max && Pop[ii].local_search == true)
                                                        {
                                                            max = Pop[ii].F;
                                                            max_index = ii;
                                                        }
                                                    }
                                                    int randx = Rand.Next(0, NPop);
                                                    Pop[randx].x = (double[])Best_nodes.x.Clone();
                                                    Pop[randx].F = Best_nodes.F;
                                                }
                                                LS Best_Indi_BLS = new LS();
                                                Best_Indi_BLS.F = Pop[max_index].F;
                                                Best_Indi_BLS.x = (double[])Pop[max_index].x.Clone();
                                                Best_Indi_BLS.visited = Pop[max_index].visited;
                                                Best_Indi_BLS.rho = Pop[max_index].rho;
                                                if (Pop[max_index].visited == 1)
                                                    Best_Indi_BLS.bias = (double[])Pop[max_index].bias.Clone();
                                                best_index = max_index;
                                                // Initialize Bias and Rho
                                                if (Best_Indi_BLS.visited == 0)
                                                {
                                                    Best_Indi_BLS.bias = new double[N];
                                                    Best_Indi_BLS.rho = 0.01;
                                                }
                                                // Local Search Operator
                                                LS Best_Indi_ALS = new LS();
                                                Best_Indi_ALS.x = (double[])Best_Indi_BLS.x.Clone();
                                                Best_Indi_ALS.F = Best_Indi_BLS.F;
                                                Best_Indi_ALS.rho = Best_Indi_BLS.rho;
                                                Best_Indi_ALS.visited = Best_Indi_BLS.visited;
                                                Best_Indi_ALS.bias = (double[])Best_Indi_BLS.bias.Clone();
                                                NFC.SolisWets(ref Best_Indi_ALS.x, ref Best_Indi_ALS.F, N, ref Best_Indi_ALS.bias, ref Best_Indi_ALS.rho, I_str, Rand, OF, Perm, Orthogonal_Matrix, M, o);
                                                if (Best_Indi_BLS.F - Best_Indi_ALS.F == 0)
                                                {
                                                    Pop[best_index].local_search = false;
                                                    //  Pop[best_index].visited = 1;
                                                }
                                                else
                                                {
                                                    Pop[best_index].x = (double[])Best_Indi_ALS.x.Clone();
                                                    Pop[best_index].visited = 1;
                                                    Pop[best_index].F = Best_Indi_ALS.F;
                                                    Pop[best_index].rho = Best_Indi_ALS.rho;
                                                    Pop[best_index].bias = (double[])Best_Indi_ALS.bias.Clone();
                                                    Pop[best_index].local_search = true;
                                                }

                                                if (Best_nodes.F < Pop[best_index].F)
                                                {
                                                    Best_nodes.F = Pop[best_index].F;
                                                    Best_nodes.x = (double[])Pop[best_index].x.Clone();
                                                }

                                                if ((it % 2) == 0)
                                                    Costs[tt, kk++] = Best_nodes.F;
                                            }
                                            costs[tt] = Best_nodes.F;

                                        }


                                        sw = new StreamWriter(result);
                                        for (int tt = 0; tt < Ntest; tt++)
                                        {
                                            for (int ii = 0; ii < Iteration / 2; ii++)
                                            {
                                                sw.Write(Costs[tt, ii].ToString() + " ");
                                            }
                                            sw.WriteLine();
                                        }
                                        sw.Close();
                                    }
                                    //}
                                }
                            }

                            #endregion
                            #region 3LMA_test
                            else if (method == "3LMA_test")
                            {
                                O = NPop * Convert.ToInt32(NI.Text);
                                double[] N_LS_GL = { 0.3, 0.1, 0.3, 0.9, 0.5, 0.5, 0.7, 0.3, 0.9, 0.1, 0.3, 0.9, 0.7, 0, 0.9, 0.3, 0.7, 0.5, 0.9, 0.3, 0.9, 0.3 };
                                double[] Mu_rate_s = { 0.1, 0.8, 0.5, 0.8, 0.8, 1, 1, 1, 0.8, 1, 0.8, 0.8, 0.8, 0, 1, 1, 0.8, 0.8, 1, 0.5, 1, 0.01 };
                                //double[] N_LS_GL = { 0.5, 0.7 };
                                //double[] Mu_rate_s = { 0.01, 0.1, 0.5 };
                                // double[] Jrs = { 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };
                                //string[] Structs = { "Ring", "Star", "Btree", "Bipartie", "Cluster", "Grid", "Ladder", "Cross-Ladder", "Cellular" };
                                string[] Structs = { "Star" };
                                //for (int Portion_index = 0; Portion_index < 5; Portion_index++)
                                //    for (int M_index = 0; M_index < 1; M_index++)
                                //    {

                                bool OBL = true;
                                //for (int OB = 0; OB < 10; OB++)
                                //{
                                for (int SS = 0; SS < 1; SS++)
                                {
                                    string Struct = Structs[SS];
                                    //double[] N_LS_GL = { 0.3, 0.5, 0.9, 0.5, 0.5, 0.3, 0.3, 0.5, 0.5, 0.3, 0.7, 0.5, 0.5, 0, 0.5, 0.7, 0.9, 0.9, 0.9, 0.7, 0.5, 0.7 };
                                    //double[] Mu_rate_s = { 0.1, 0.00001, 5, 0.1, 0.1, 0.5, 0.01, 0.1, 5, 0.5, 0.001, 5, 0.1, 0, 0.1, 5, 0.01, 0.5, 0.01, 0.1, 0.5, 0.5 };
                                    double LS_LG = N_LS_GL[OF - 1];
                                    double Mu_rate = Mu_rate_s[OF - 1];
                                    Iteration = Convert.ToInt32((O * (1 - LS_LG)) / NPop) + 1;
                                    int I_str = Convert.ToInt32(O * LS_LG) / Iteration;
                                    int best_index = 0;
                                    bool local_search = false;
                                    double Jr = 0.3;
                                    double[,] Costs = new double[Ntest, Iteration / 2 + 1];
                                    //ss2 = "fit_N_" + N.ToString() + "_Mr_" + Mu_rate.ToString() + "LSLG" + LS_LG.ToString() + "_J" + Jr.ToString() + "_Pop_" + NPop.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "Struct_" + Struct + "_OF" + OF.ToString() + ".txt";
                                    if (!File.Exists(result))
                                    {
                                        for (int tt = 0; tt < Ntest; tt++)
                                        {

                                            //counter_JR = 1;
                                            LS[] Pop = new LS[NPop];
                                            LS[] OPop = new LS[NPop];
                                            NF Best_nodes = new NF();
                                            NF[] Breed = new NF[NPop];
                                            Best_nodes.F = -1E250;
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                Pop[ii].x = NFC.node_constructing(N, Rand);
                                                if (OBL)
                                                    OPop[ii].x = NFC.Opposition_Learning(Pop[ii].x);
                                                if (OF < 16 || OF == 22)
                                                    Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N);
                                                else
                                                    Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                if (OBL)
                                                {
                                                    if (OF < 16 || OF == 22)
                                                        OPop[ii].F = NFC.compute_fitness(OPop[ii].x, OF, N);
                                                    else
                                                        OPop[ii].F = NFC.compute_fitness(OPop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);

                                                    if (OPop[ii].F > Pop[ii].F)
                                                    {
                                                        Pop[ii].F = OPop[ii].F;
                                                        Pop[ii].x = (double[])OPop[ii].x.Clone();
                                                    }
                                                }
                                                Pop[ii].local_search = true;

                                            }
                                            int OBL_iter = 0;
                                            int kk = 0;
                                            for (int it = 0; it < Iteration; it++)
                                            {
                                                //perform Steady State GA
                                                #region generate Breeds
                                                double max = -1E250; int max_index = 0;
                                                local_search = false;
                                                int imax = -1; int imin = -1; double max_b = -1E250; double min_p = 1E250;
                                                for (int ii = 0; ii < NPop; ii++)
                                                {
                                                    int parent1 = Rand.Next(NPop);
                                                    int[] S = NFC.Structures(ii, NPop, Struct);
                                                    int index2 = Rand.Next(S.Length);
                                                    //int parent2 = S[index2];
                                                    int parent2 = NFC.negative_assortative_mating_strategy(Pop[parent1].x, Pop, Rand, S);
                                                    Breed[ii].x = NFC.BLX_Crossover(Pop[parent1].x, Pop[parent2].x);
                                                    if (Rand.NextDouble() < Mu_rate)
                                                        NFC.Mutation_BGA2(ref Breed[ii], Pop, N, Rand);
                                                    if (OF < 16 || OF == 22)
                                                        Breed[ii].F = NFC.compute_fitness(Breed[ii].x, OF, N);
                                                    else
                                                        Breed[ii].F = NFC.compute_fitness(Breed[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                    Breed[ii].x = NFC.Node_treatment(Breed[ii].x, N, Rand);
                                                    ////////////////////////////////////////////
                                                    if (Breed[ii].F > max_b)
                                                    {
                                                        imax = ii;
                                                        max_b = Breed[ii].F;
                                                    }
                                                    if (Pop[ii].F < min_p)
                                                    {
                                                        imin = ii;
                                                        min_p = Pop[ii].F;
                                                    }
                                                }
                                                if (max_b > min_p)
                                                {
                                                    Pop[imin].F = Breed[imax].F;
                                                    Pop[imin].x = (double[])Breed[imax].x.Clone();
                                                }
                                                for (int ii = 0; ii < NPop; ii++)
                                                {
                                                    if (OBL)
                                                    {
                                                        if (Jr > Rand.NextDouble())
                                                        {

                                                            OPop[ii].x = NFC.generating_jumping(Pop, Pop[ii].x, NPop, N, ii);
                                                            if (OF < 16 || OF == 22)
                                                                OPop[ii].F = NFC.compute_fitness(OPop[ii].x, OF, N);
                                                            else
                                                                OPop[ii].F = NFC.compute_fitness(OPop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                            if (Pop[ii].F < OPop[ii].F)
                                                            {
                                                                Pop[ii].x = (double[])OPop[ii].x.Clone();
                                                                Pop[ii].F = OPop[ii].F;
                                                            }

                                                            OBL_iter++;
                                                            if (OBL_iter == NPop)
                                                            {
                                                                OBL_iter = 0;
                                                                it++;
                                                            }
                                                        }
                                                    }
                                                    if (Pop[ii].F > max && Pop[ii].local_search == true)
                                                    {
                                                        max = Pop[ii].F;
                                                        max_index = ii;
                                                        local_search = true;
                                                    }
                                                }
                                                #endregion
                                                // Replacement Strategy

                                                //
                                                //Pick the best individual of Population
                                                if (local_search == false)
                                                {

                                                    for (int ii = 0; ii < NPop; ii++)
                                                    {
                                                        Pop[ii].x = NFC.node_constructing(N, Rand);
                                                        if (OF < 16 || OF == 22)
                                                            Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N);
                                                        else
                                                            Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                        Pop[ii].x = NFC.Node_treatment(Pop[ii].x, N, Rand);
                                                        Pop[ii].local_search = true;
                                                        Pop[ii].visited = 0;
                                                        if (Pop[ii].F > max && Pop[ii].local_search == true)
                                                        {
                                                            max = Pop[ii].F;
                                                            max_index = ii;
                                                        }
                                                    }
                                                    int randx = Rand.Next(0, NPop);
                                                    Pop[randx].x = (double[])Best_nodes.x.Clone();
                                                    Pop[randx].F = Best_nodes.F;
                                                }
                                                LS Best_Indi_BLS = new LS();
                                                Best_Indi_BLS.F = Pop[max_index].F;
                                                Best_Indi_BLS.x = (double[])Pop[max_index].x.Clone();
                                                Best_Indi_BLS.visited = Pop[max_index].visited;
                                                Best_Indi_BLS.rho = Pop[max_index].rho;
                                                if (Pop[max_index].visited == 1)
                                                    Best_Indi_BLS.bias = (double[])Pop[max_index].bias.Clone();
                                                best_index = max_index;
                                                // Initialize Bias and Rho
                                                if (Best_Indi_BLS.visited == 0)
                                                {
                                                    Best_Indi_BLS.bias = new double[N];
                                                    Best_Indi_BLS.rho = 0.01;
                                                }
                                                // Local Search Operator
                                                LS Best_Indi_ALS = new LS();
                                                Best_Indi_ALS.x = (double[])Best_Indi_BLS.x.Clone();
                                                Best_Indi_ALS.F = Best_Indi_BLS.F;
                                                Best_Indi_ALS.rho = Best_Indi_BLS.rho;
                                                Best_Indi_ALS.visited = Best_Indi_BLS.visited;
                                                Best_Indi_ALS.bias = (double[])Best_Indi_BLS.bias.Clone();
                                                NFC.SolisWets(ref Best_Indi_ALS.x, ref Best_Indi_ALS.F, N, ref Best_Indi_ALS.bias, ref Best_Indi_ALS.rho, I_str, Rand, OF, Perm, Orthogonal_Matrix, M, o);
                                                if (Best_Indi_BLS.F - Best_Indi_ALS.F == 0)
                                                {
                                                    Pop[best_index].local_search = false;
                                                    //  Pop[best_index].visited = 1;
                                                }
                                                else
                                                {
                                                    Pop[best_index].x = (double[])Best_Indi_ALS.x.Clone();
                                                    Pop[best_index].visited = 1;
                                                    Pop[best_index].F = Best_Indi_ALS.F;
                                                    Pop[best_index].rho = Best_Indi_ALS.rho;
                                                    Pop[best_index].bias = (double[])Best_Indi_ALS.bias.Clone();
                                                    Pop[best_index].local_search = true;
                                                }

                                                if (Best_nodes.F < Pop[best_index].F)
                                                {
                                                    Best_nodes.F = Pop[best_index].F;
                                                    Best_nodes.x = (double[])Pop[best_index].x.Clone();
                                                }

                                                //if ((it % 2) == 0)
                                                //    Costs[tt, kk++] = Best_nodes.F;
                                            }
                                            costs[tt] = Best_nodes.F;

                                        }


                                        sw = new StreamWriter(result);
                                        for (int tt = 0; tt < Ntest; tt++)
                                        {
                                            //for (int ii = 0; ii < Iteration / 2; ii++)
                                            //{
                                            //if(Costs[tt, ii]!=0)
                                            sw.Write(costs[tt].ToString() + " ");
                                            //}
                                            //sw.WriteLine();
                                        }
                                        sw.Close();
                                    }
                                    //}
                                }
                            }

                            #endregion
                            #region 3LMAOriginal
                            else if (method == "2LMA")
                            {
                                O = NPop * Convert.ToInt32(NI.Text);
                                double[] N_LS_GL = { 0.3, 0.1, 0.3, 0.9, 0.5, 0.5, 0.7, 0.3, 0.9, 0.1, 0.3, 0.9, 0.7, 0, 0.9, 0.3, 0.7, 0.5, 0.9, 0.3, 0.9, 0.3 };
                                double[] Mu_rate_s = { 0.1, 0.8, 0.5, 0.8, 0.8, 1, 1, 1, 0.8, 1, 0.8, 0.8, 0.8, 0, 1, 1, 0.8, 0.8, 1, 0.5, 1, 0.01 };
                                //double[] N_LS_GL = { 0.5, 0.7 };
                                //double[] Mu_rate_s = { 0.01, 0.1, 0.5 };
                                // double[] Jrs = { 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };
                                //string[] Structs = { "Ring", "Star", "Btree", "Bipartie", "Cluster", "Grid", "Ladder", "Cross-Ladder", "Cellular" };
                                string[] Structs = { "Star" };
                                //for (int Portion_index = 0; Portion_index < 5; Portion_index++)
                                //    for (int M_index = 0; M_index < 1; M_index++)
                                //    {

                                bool OBL = true;
                                //for (int OB = 0; OB < 10; OB++)
                                //{
                                for (int SS = 0; SS < 1; SS++)
                                {
                                    string Struct = Structs[SS];
                                    //double[] N_LS_GL = { 0.3, 0.5, 0.9, 0.5, 0.5, 0.3, 0.3, 0.5, 0.5, 0.3, 0.7, 0.5, 0.5, 0, 0.5, 0.7, 0.9, 0.9, 0.9, 0.7, 0.5, 0.7 };
                                    //double[] Mu_rate_s = { 0.1, 0.00001, 5, 0.1, 0.1, 0.5, 0.01, 0.1, 5, 0.5, 0.001, 5, 0.1, 0, 0.1, 5, 0.01, 0.5, 0.01, 0.1, 0.5, 0.5 };
                                    double LS_LG = N_LS_GL[OF - 1];
                                    double Mu_rate = Mu_rate_s[OF - 1];
                                    Iteration = Convert.ToInt32((O * (1 - LS_LG)) / NPop) + 1;
                                    int I_str = Convert.ToInt32(O * LS_LG) / Iteration;
                                    int best_index = 0;
                                    bool local_search = false;
                                    double Jr = 0.3;
                                    double[,] Costs = new double[Ntest, Iteration / 2 + 1];
                                    //ss2 = "fit_N_" + N.ToString() + "_Mr_" + Mu_rate.ToString() + "LSLG" + LS_LG.ToString() + "_J" + Jr.ToString() + "_Pop_" + NPop.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "Struct_" + Struct + "_OF" + OF.ToString() + ".txt";
                                    if (!File.Exists(result))
                                    {
                                        for (int tt = 0; tt < Ntest; tt++)
                                        {

                                            //counter_JR = 1;
                                            LS[] Pop = new LS[NPop];
                                            LS[] OPop = new LS[NPop];
                                            NF Best_nodes = new NF();
                                            NF[] Breed = new NF[NPop];
                                            Best_nodes.F = -1E250;
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                Pop[ii].x = NFC.node_constructing(N, Rand);
                                                if (OBL)
                                                    OPop[ii].x = NFC.Opposition_Learning(Pop[ii].x);
                                                if (OF < 16 || OF == 22)
                                                    Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N);
                                                else
                                                    Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                if (OBL)
                                                {
                                                    if (OF < 16 || OF == 22)
                                                        OPop[ii].F = NFC.compute_fitness(OPop[ii].x, OF, N);
                                                    else
                                                        OPop[ii].F = NFC.compute_fitness(OPop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);

                                                    if (OPop[ii].F > Pop[ii].F)
                                                    {
                                                        Pop[ii].F = OPop[ii].F;
                                                        Pop[ii].x = (double[])OPop[ii].x.Clone();
                                                    }
                                                }
                                                Pop[ii].local_search = true;

                                            }

                                            int kk = 0;
                                            int OBL_iter = 0;
                                            for (int it = 0; it < Iteration; it++)
                                            {
                                                //perform Steady State GA
                                                #region generate Breeds
                                                double max = -1E250; int max_index = 0;
                                                local_search = false;
                                                for (int ii = 0; ii < NPop; ii++)
                                                {
                                                    int parent1 = Rand.Next(NPop);
                                                    int[] S = NFC.Structures(ii, NPop, Struct);
                                                    int index2 = Rand.Next(S.Length);
                                                    //int parent2 = S[index2];
                                                    int parent2 = NFC.negative_assortative_mating_strategy(Pop[parent1].x, Pop, Rand, S);
                                                    Breed[ii].x = NFC.BLX_Crossover(Pop[parent1].x, Pop[parent2].x);
                                                    if (Rand.NextDouble() < Mu_rate)
                                                        NFC.Mutation_BGA2(ref Breed[ii], Pop, N, Rand);
                                                    if (OF < 16 || OF == 22)
                                                        Breed[ii].F = NFC.compute_fitness(Breed[ii].x, OF, N);
                                                    else
                                                        Breed[ii].F = NFC.compute_fitness(Breed[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                    Breed[ii].x = NFC.Node_treatment(Breed[ii].x, N, Rand);
                                                    ////////////////////////////////////////////
                                                    if (Breed[ii].F > Pop[ii].F)
                                                    {
                                                        Pop[ii].x = (double[])Breed[ii].x.Clone();
                                                        Pop[ii].F = Breed[ii].F;
                                                    }
                                                }
                                                for (int ii = 0; ii < NPop; ii++)
                                                {
                                                    if (OBL)
                                                    {
                                                        if (Jr > Rand.NextDouble())
                                                        {

                                                            OPop[ii].x = NFC.generating_jumping(Pop, Pop[ii].x, NPop, N, ii);
                                                            if (OF < 16 || OF == 22)
                                                                OPop[ii].F = NFC.compute_fitness(OPop[ii].x, OF, N);
                                                            else
                                                                OPop[ii].F = NFC.compute_fitness(OPop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                            if (Pop[ii].F < OPop[ii].F)
                                                            {
                                                                Pop[ii].x = (double[])OPop[ii].x.Clone();
                                                                Pop[ii].F = OPop[ii].F;
                                                            }
                                                            OBL_iter++;
                                                            if (OBL_iter == NPop)
                                                            {
                                                                OBL_iter = 0;
                                                                it++;
                                                            }
                                                        }
                                                    }
                                                    if (Pop[ii].F > max && Pop[ii].local_search == true)
                                                    {
                                                        max = Pop[ii].F;
                                                        max_index = ii;
                                                        local_search = true;
                                                    }
                                                }
                                                #endregion
                                                // Replacement Strategy

                                                //
                                                //Pick the best individual of Population
                                                if (local_search == false)
                                                {

                                                    for (int ii = 0; ii < NPop; ii++)
                                                    {
                                                        Pop[ii].x = NFC.node_constructing(N, Rand);
                                                        if (OF < 16 || OF == 22)
                                                            Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N);
                                                        else
                                                            Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                        Pop[ii].x = NFC.Node_treatment(Pop[ii].x, N, Rand);
                                                        Pop[ii].local_search = true;
                                                        Pop[ii].visited = 0;
                                                        if (Pop[ii].F > max && Pop[ii].local_search == true)
                                                        {
                                                            max = Pop[ii].F;
                                                            max_index = ii;
                                                        }
                                                    }
                                                    int randx = Rand.Next(0, NPop);
                                                    Pop[randx].x = (double[])Best_nodes.x.Clone();
                                                    Pop[randx].F = Best_nodes.F;
                                                }
                                                LS Best_Indi_BLS = new LS();
                                                Best_Indi_BLS.F = Pop[max_index].F;
                                                Best_Indi_BLS.x = (double[])Pop[max_index].x.Clone();
                                                Best_Indi_BLS.visited = Pop[max_index].visited;
                                                Best_Indi_BLS.rho = Pop[max_index].rho;
                                                if (Pop[max_index].visited == 1)
                                                    Best_Indi_BLS.bias = (double[])Pop[max_index].bias.Clone();
                                                best_index = max_index;
                                                // Initialize Bias and Rho
                                                if (Best_Indi_BLS.visited == 0)
                                                {
                                                    Best_Indi_BLS.bias = new double[N];
                                                    Best_Indi_BLS.rho = 0.01;
                                                }
                                                // Local Search Operator
                                                LS Best_Indi_ALS = new LS();
                                                Best_Indi_ALS.x = (double[])Best_Indi_BLS.x.Clone();
                                                Best_Indi_ALS.F = Best_Indi_BLS.F;
                                                Best_Indi_ALS.rho = Best_Indi_BLS.rho;
                                                Best_Indi_ALS.visited = Best_Indi_BLS.visited;
                                                Best_Indi_ALS.bias = (double[])Best_Indi_BLS.bias.Clone();
                                                NFC.SolisWets(ref Best_Indi_ALS.x, ref Best_Indi_ALS.F, N, ref Best_Indi_ALS.bias, ref Best_Indi_ALS.rho, I_str, Rand, OF, Perm, Orthogonal_Matrix, M, o);
                                                if (Best_Indi_BLS.F - Best_Indi_ALS.F == 0)
                                                {
                                                    Pop[best_index].local_search = false;
                                                    //  Pop[best_index].visited = 1;
                                                }
                                                else
                                                {
                                                    Pop[best_index].x = (double[])Best_Indi_ALS.x.Clone();
                                                    Pop[best_index].visited = 1;
                                                    Pop[best_index].F = Best_Indi_ALS.F;
                                                    Pop[best_index].rho = Best_Indi_ALS.rho;
                                                    Pop[best_index].bias = (double[])Best_Indi_ALS.bias.Clone();
                                                    Pop[best_index].local_search = true;
                                                }

                                                if (Best_nodes.F < Pop[best_index].F)
                                                {
                                                    Best_nodes.F = Pop[best_index].F;
                                                    Best_nodes.x = (double[])Pop[best_index].x.Clone();
                                                }

                                                //if ((it % 2) == 0)
                                                //    Costs[tt, kk++] = Best_nodes.F;
                                            }
                                            costs[tt] = Best_nodes.F;

                                        }


                                        sw = new StreamWriter(result);
                                        for (int tt = 0; tt < Ntest; tt++)
                                        {
                                            //    for (int ii = 0; ii < Iteration / 2; ii++)
                                            //    {
                                            sw.Write(costs[tt].ToString() + " ");
                                            //}
                                            //    sw.WriteLine();
                                        }
                                        sw.Close();
                                    }
                                    //}
                                }
                            }
                            #endregion
                            #region 12LMA
                            else if (method == "12LMA")
                            {
                                O = NPop * Convert.ToInt32(NI.Text);
                                //double[] N_LS_GL = { 0.3, 0.1, 0.3, 0.9, 0.5, 0.5, 0.7, 0.3, 0.9, 0.1, 0.3, 0.9, 0.7, 0, 0.9, 0.3, 0.7, 0.5, 0.9, 0.3, 0.9, 0.3 };
                                double[] Mu_rate_s = { 0.1, 0.8, 0.5, 0.8, 0.8, 1, 1, 1, 0.8, 1, 0.8, 0.8, 0.8, 0, 1, 1, 0.8, 0.8, 1, 0.5, 1, 0.01 };
                                //double[] N_LS_GL = { 0.5, 0.7 };
                                //double[] Mu_rate_s = { 0.01, 0.1, 0.5 };
                                // double[] Jrs = { 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };
                                //string[] Structs = { "Ring", "Star", "Btree", "Bipartie", "Cluster", "Grid", "Ladder", "Cross-Ladder", "Cellular" };
                                string[] Structs = { "Star" };
                                //for (int Portion_index = 0; Portion_index < 5; Portion_index++)
                                //    for (int M_index = 0; M_index < 1; M_index++)
                                //    {

                                bool OBL = true;
                                //for (int OB = 0; OB < 10; OB++)
                                //{
                                for (int SS = 0; SS < 1; SS++)
                                {
                                    string Struct = Structs[SS];
                                    //double[] N_LS_GL = { 0.3, 0.5, 0.9, 0.5, 0.5, 0.3, 0.3, 0.5, 0.5, 0.3, 0.7, 0.5, 0.5, 0, 0.5, 0.7, 0.9, 0.9, 0.9, 0.7, 0.5, 0.7 };
                                    //double[] Mu_rate_s = { 0.1, 0.00001, 5, 0.1, 0.1, 0.5, 0.01, 0.1, 5, 0.5, 0.001, 5, 0.1, 0, 0.1, 5, 0.01, 0.5, 0.01, 0.1, 0.5, 0.5 };
                                    // double LS_LG = N_LS_GL[OF - 1];
                                    double Mu_rate = Mu_rate_s[OF - 1];
                                    // Iteration = Convert.ToInt32((O * (1 - LS_LG)) / NPop) + 1;
                                    // int I_str = Convert.ToInt32(O * LS_LG) / Iteration;
                                    int best_index = 0;
                                    double Jr = 0.3;
                                    double max = -1E290;
                                    double[,] Costs = new double[Ntest, Iteration / 2 + 1];
                                    //ss2 = "fit_N_" + N.ToString() + "_Mr_" + Mu_rate.ToString() + "_J" + Jr.ToString() + "_Pop_" + NPop.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "Struct_" + Struct + "_OF" + OF.ToString() + ".txt";
                                    if (!File.Exists(result))
                                    {
                                        for (int tt = 0; tt < Ntest; tt++)
                                        {

                                            //counter_JR = 1;
                                            LS[] Pop = new LS[NPop];
                                            LS[] OPop = new LS[NPop];
                                            NF Best_nodes = new NF();
                                            NF[] Breed = new NF[NPop];
                                            Best_nodes.F = -1E250;
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                Pop[ii].x = NFC.node_constructing(N, Rand);
                                                if (OBL)
                                                    OPop[ii].x = NFC.Opposition_Learning(Pop[ii].x);
                                                if (OF < 16 || OF == 22)
                                                    Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N);
                                                else
                                                    Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                if (OBL)
                                                {
                                                    if (OF < 16 || OF == 22)
                                                        OPop[ii].F = NFC.compute_fitness(OPop[ii].x, OF, N);
                                                    else
                                                        OPop[ii].F = NFC.compute_fitness(OPop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);

                                                    if (OPop[ii].F > Pop[ii].F)
                                                    {
                                                        Pop[ii].F = OPop[ii].F;
                                                        Pop[ii].x = (double[])OPop[ii].x.Clone();
                                                    }
                                                }
                                                Pop[ii].local_search = true;

                                            }

                                            int kk = 0;
                                            int OBL_iter = 0;
                                            for (int it = 0; it < Iteration; it++)
                                            {
                                                //perform Steady State GA
                                                #region generate Breeds
                                                int imax = -1; int imin = -1; double max_b = -1E250; double min_p = 1E250;
                                                for (int ii = 0; ii < NPop; ii++)
                                                {
                                                    int parent1 = Rand.Next(NPop);
                                                    int[] S = NFC.Structures(ii, NPop, Struct);
                                                    int index2 = Rand.Next(S.Length);
                                                    //int parent2 = S[index2];
                                                    int parent2 = NFC.negative_assortative_mating_strategy(Pop[parent1].x, Pop, Rand, S);
                                                    Breed[ii].x = NFC.BLX_Crossover(Pop[parent1].x, Pop[parent2].x);
                                                    if (Rand.NextDouble() < Mu_rate)
                                                        NFC.Mutation_BGA2(ref Breed[ii], Pop, N, Rand);
                                                    if (OF < 16 || OF == 22)
                                                        Breed[ii].F = NFC.compute_fitness(Breed[ii].x, OF, N);
                                                    else
                                                        Breed[ii].F = NFC.compute_fitness(Breed[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                    Breed[ii].x = NFC.Node_treatment(Breed[ii].x, N, Rand);
                                                    ////////////////////////////////////////////
                                                    if (Breed[ii].F > max_b)
                                                    {
                                                        imax = ii;
                                                        max_b = Breed[ii].F;
                                                    }
                                                    if (Pop[ii].F < min_p)
                                                    {
                                                        imin = ii;
                                                        min_p = Pop[ii].F;
                                                    }
                                                }
                                                if (max_b > min_p)
                                                {
                                                    Pop[imin].F = Breed[imax].F;
                                                    Pop[imin].x = (double[])Breed[imax].x.Clone();
                                                }
                                                for (int ii = 0; ii < NPop; ii++)
                                                {
                                                    if (OBL)
                                                    {
                                                        if (Jr > Rand.NextDouble())
                                                        {

                                                            OPop[ii].x = NFC.generating_jumping(Pop, Pop[ii].x, NPop, N, ii);
                                                            if (OF < 16 || OF == 22)
                                                                OPop[ii].F = NFC.compute_fitness(OPop[ii].x, OF, N);
                                                            else
                                                                OPop[ii].F = NFC.compute_fitness(OPop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                            if (Pop[ii].F < OPop[ii].F)
                                                            {
                                                                Pop[ii].x = (double[])OPop[ii].x.Clone();
                                                                Pop[ii].F = OPop[ii].F;
                                                            }
                                                            OBL_iter++;
                                                            if (OBL_iter == NPop)
                                                            {
                                                                OBL_iter = 0;
                                                                it++;
                                                            }
                                                        }
                                                    }
                                                    if (Pop[ii].F > max)
                                                    {
                                                        max = Pop[ii].F;
                                                    }
                                                }
                                                #endregion

                                                if (Best_nodes.F < Pop[best_index].F)
                                                {
                                                    Best_nodes.F = Pop[best_index].F;
                                                    Best_nodes.x = (double[])Pop[best_index].x.Clone();
                                                }

                                                //if ((it % 2) == 0)
                                                //    Costs[tt, kk++] = Best_nodes.F;
                                            }
                                            costs[tt] = Best_nodes.F;

                                        }


                                        sw = new StreamWriter(result);
                                        for (int tt = 0; tt < Ntest; tt++)
                                        {
                                            //for (int ii = 0; ii < Iteration / 2; ii++)
                                            //{
                                            //if (Costs[tt, ii] != 0)
                                            sw.Write(costs[tt].ToString() + " ");
                                            //}
                                            //sw.WriteLine();
                                        }
                                        sw.Close();
                                    }
                                    //}
                                }
                            }

                            #endregion
                            #region 1LMA
                            else if (method == "1LMA")
                            {
                                O = NPop * Convert.ToInt32(NI.Text);
                                //double[] N_LS_GL = { 0.3, 0.1, 0.3, 0.9, 0.5, 0.5, 0.7, 0.3, 0.9, 0.1, 0.3, 0.9, 0.7, 0, 0.9, 0.3, 0.7, 0.5, 0.9, 0.3, 0.9, 0.3 };
                                double[] Mu_rate_s = { 0.1, 0.8, 0.5, 0.8, 0.8, 1, 1, 1, 0.8, 1, 0.8, 0.8, 0.8, 0, 1, 1, 0.8, 0.8, 1, 0.5, 1, 0.01 };
                                string[] Structs = { "Star" };
                                bool OBL = true;
                                for (int SS = 0; SS < 1; SS++)
                                {
                                    string Struct = Structs[SS];
                                    //double[] N_LS_GL = { 0.3, 0.5, 0.9, 0.5, 0.5, 0.3, 0.3, 0.5, 0.5, 0.3, 0.7, 0.5, 0.5, 0, 0.5, 0.7, 0.9, 0.9, 0.9, 0.7, 0.5, 0.7 };
                                    //double[] Mu_rate_s = { 0.1, 0.00001, 5, 0.1, 0.1, 0.5, 0.01, 0.1, 5, 0.5, 0.001, 5, 0.1, 0, 0.1, 5, 0.01, 0.5, 0.01, 0.1, 0.5, 0.5 };
                                    // double LS_LG = N_LS_GL[OF - 1];
                                    double Mu_rate = Mu_rate_s[OF - 1];
                                    // Iteration = Convert.ToInt32((O * (1 - LS_LG)) / NPop) + 1;
                                    // int I_str = Convert.ToInt32(O * LS_LG) / Iteration;
                                    int best_index = 0;
                                    double Jr = 0.3;
                                    double max = -1E290;
                                    double[,] Costs = new double[Ntest, Iteration / 2 + 1];
                                    //ss2 = "fit_N_" + N.ToString() + "_Mr_" + Mu_rate.ToString() + "_J" + Jr.ToString() + "_Pop_" + NPop.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "Struct_" + Struct + "_OF" + OF.ToString() + ".txt";
                                    if (!File.Exists(result))
                                    {
                                        for (int tt = 0; tt < Ntest; tt++)
                                        {

                                            //counter_JR = 1;
                                            LS[] Pop = new LS[NPop];
                                            LS[] OPop = new LS[NPop];
                                            NF Best_nodes = new NF();
                                            NF[] Breed = new NF[NPop];
                                            Best_nodes.F = -1E250;
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                Pop[ii].x = NFC.node_constructing(N, Rand);
                                                OPop[ii].x = NFC.Opposition_Learning(Pop[ii].x);
                                                if (OF < 16 || OF == 22)
                                                    Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N);
                                                else
                                                    Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);


                                            }

                                            int kk = 0;
                                            for (int it = 0; it < Iteration; it++)
                                            {
                                                //perform Steady State GA
                                                #region generate Breeds
                                                int imax = -1; int imin = -1; double max_b = -1E250; double min_p = 1E250;
                                                for (int ii = 0; ii < NPop; ii++)
                                                {
                                                    if (Best_nodes.F < Pop[ii].F)
                                                    {
                                                        Best_nodes.F = Pop[best_index].F;
                                                        Best_nodes.x = (double[])Pop[best_index].x.Clone();
                                                    }
                                                    int parent1 = Rand.Next(NPop);
                                                    int[] S = NFC.Structures(ii, NPop, Struct);
                                                    int index2 = Rand.Next(S.Length);
                                                    //int parent2 = S[index2];
                                                    int parent2 = NFC.negative_assortative_mating_strategy(Pop[parent1].x, Pop, Rand, S);
                                                    Breed[ii].x = NFC.BLX_Crossover(Pop[parent1].x, Pop[parent2].x);
                                                    if (Rand.NextDouble() < Mu_rate)
                                                        NFC.Mutation_BGA2(ref Breed[ii], Pop, N, Rand);
                                                    if (OF < 16 || OF == 22)
                                                        Breed[ii].F = NFC.compute_fitness(Breed[ii].x, OF, N);
                                                    else
                                                        Breed[ii].F = NFC.compute_fitness(Breed[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                    Breed[ii].x = NFC.Node_treatment(Breed[ii].x, N, Rand);
                                                    ////////////////////////////////////////////
                                                    if (Breed[ii].F > max_b)
                                                    {
                                                        imax = ii;
                                                        max_b = Breed[ii].F;
                                                    }
                                                    if (Pop[ii].F < min_p)
                                                    {
                                                        imin = ii;
                                                        min_p = Pop[ii].F;
                                                    }
                                                }
                                                if (max_b > min_p)
                                                {
                                                    Pop[imin].F = Breed[imax].F;
                                                    Pop[imin].x = (double[])Breed[imax].x.Clone();
                                                }
                                                #endregion

                                                if ((it % 2) == 0)
                                                    Costs[tt, kk++] = Best_nodes.F;
                                            }
                                            costs[tt] = Best_nodes.F;

                                        }


                                        sw = new StreamWriter(result);
                                        for (int tt = 0; tt < Ntest; tt++)
                                        {
                                            for (int ii = 0; ii < Iteration / 2; ii++)
                                            {
                                                sw.Write(Costs[tt, ii].ToString() + " ");
                                            }
                                            sw.WriteLine();
                                        }
                                        sw.Close();
                                    }
                                    //}
                                }
                            }

                            #endregion
                            #region 3SOME
                            if (method == "3SOME")
                            {

                                int k = 4;
                                //double[]inheritance_factors={0.01,0.1,0.3,0.6,0.9};
                                //double[] deltas = { 0.001, 0.01, 0.1, 0.4, 0.6 };
                                //double[] inheritance_factors = { 0.6, 0.6, 0.9, 0.6, 0.3, 0.9 };
                                //double[] deltas = { 0.01, 0.4, 0.001, 0.4, 0.001, 0.001 };
                                double[] inheritance_factors = { 0.3, 0.6, 0.3, 0.6, 0.6, 0.9, 0.3, 0.6, 0.6, 0.3, 0.6, 0.9, 0.9, 0, 0.6, 0.9, 0.6, 0.1, 0.3, 0.9, 0.01, 0.6 };
                                double[] deltas = { 0.001, 0.01, 0.4, 0.001, 0.001, 0.01, 0.4, 0.01, 0.001, 0.1, 0.01, 0.01, 0.001, 0, 0.001, 0.001, 0.01, 0.01, 0.001, 0.1, 0.01, 0.01 };
                                //double[] rhos = { 0.01, 0.1, 0.2, 0.4, 0.6 };
                                bool updated = false;
                                //for (int inh = 0; inh < 5; inh++)
                                //    for (int del = 0; del < 5; del++)
                                //    {
                                double rho = 0.4;
                                double inheritance_factor = inheritance_factors[OF];
                                double delta = deltas[OF];
                                double root = N * inheritance_factor;
                                double Cr = 1.0 / (Math.Pow(2, 1.0 / root));
                                root = N * (1 - inheritance_factor);
                                double Crprim = 1.0 / (Math.Pow(2, 1.0 / root));
                                int local_budget = 150;
                                //ss2 = "Results\\fit_N_" + N.ToString() + "_bound_" + boundary_method.ToString() + "pop_" + NPop.ToString() + method + "_inh_" + inheritance_factor.ToString() + "_del_" + delta.ToString() + "_rho_" + rho.ToString() + "_test_" + Ntest.ToString() + "OF_" + OF.ToString() + ".txt";
                                #region file_existency
                                if (!File.Exists(result))
                                {
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        rho = 0.4;
                                        NF Elite = new NF();
                                        NF Indi = new NF();
                                        ////////////////////////////////////
                                        //////////////////////////////
                                        Elite.x = NFC.node_constructing(N, Rand);
                                        if (OF < 16 || OF == 22)
                                            Elite.F = NFC.compute_fitness(Elite.x, OF, N);
                                        else
                                            Elite.F = NFC.compute_fitness(Elite.x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                        for (int it = 0; it < Iteration * NPop; )
                                        {
                                            ///////////////////Long distance/////////////////////////////////
                                            updated = false;
                                            while (!updated && it < Iteration * NPop)
                                            {
                                                #region Long_distance
                                                int checker = 0;
                                                Indi.x = NFC.node_constructing(N, Rand);
                                                int i = Rand.Next(N);
                                                Indi.x[i] = Elite.x[i];
                                                while (Rand.NextDouble() <= Cr)
                                                {
                                                    Indi.x[i] = Elite.x[i];
                                                    i++;
                                                    if (i == N)
                                                        i = 0;
                                                    if (checker >= N)
                                                        break;
                                                    checker++;
                                                }
                                                if (OF < 16 || OF == 22)
                                                    Indi.F = NFC.compute_fitness(Indi.x, OF, N);
                                                else
                                                    Indi.F = NFC.compute_fitness(Indi.x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                if (Indi.F >= Elite.F)
                                                {
                                                    Elite.x = (double[])Indi.x.Clone();
                                                    Elite.F = Indi.F;
                                                    updated = true;
                                                }
                                                //////////////////////Function Evaluation increasing///////////////////
                                                it++;
                                                #endregion
                                            }
                                            ///////////////Middle Distance////////////////////////////
                                            while (updated && it < Iteration * NPop)
                                            {
                                                updated = false;
                                                #region Middle_distance
                                                for (int I = 0; I < k * N; I++)
                                                {
                                                    Indi.x = NFC.Initialization_Hypercube(Elite.x, N, Rand, delta);
                                                    int i = Rand.Next(N);
                                                    Indi.x[i] = Elite.x[i];
                                                    int checker = 0;
                                                    while (Rand.NextDouble() <= Crprim)
                                                    {
                                                        Indi.x[i] = Elite.x[i];
                                                        i++;
                                                        if (i == N)
                                                            i = 0;
                                                        if (checker >= N)
                                                            break;
                                                        checker++;
                                                    }
                                                    //Indi.x = NFC.Node_treatment(Indi.x, N, Rand);
                                                    Solution_Treatment(boundary_method, ref Indi.x, Rand);
                                                    if (OF < 16 || OF == 22)
                                                        Indi.F = NFC.compute_fitness(Indi.x, OF, N);
                                                    else
                                                        Indi.F = NFC.compute_fitness(Indi.x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                    if (Indi.F > Elite.F)
                                                    {
                                                        updated = true;
                                                    }
                                                    if (Indi.F >= Elite.F)
                                                    {
                                                        Elite.x = (double[])Indi.x.Clone();
                                                        Elite.F = Indi.F;

                                                    }
                                                    it++;
                                                }
                                                #endregion
                                            }
                                            /////////////////////////Short Distance////////////////////////////////////////
                                            #region Short_region
                                            NF indi2 = new NF();
                                            for (int ii = 0; ii < local_budget && it < NPop * Iteration; ii++)
                                            {
                                                Indi.x = (double[])Elite.x.Clone();
                                                indi2.x = (double[])Elite.x.Clone();
                                                for (int jj = 0; jj < N; jj++)
                                                {
                                                    if (Elite.x[jj] - rho > -1)
                                                        indi2.x[jj] = Elite.x[jj] - rho;
                                                    else
                                                        indi2.x[jj] = Elite.x[jj];
                                                    if (OF < 16 || OF == 22)
                                                        indi2.F = NFC.compute_fitness(indi2.x, OF, N);
                                                    else
                                                        indi2.F = NFC.compute_fitness(indi2.x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                    it++;
                                                    if (indi2.F >= Indi.F)
                                                    {
                                                        Indi.x = (double[])indi2.x.Clone();
                                                        Indi.F = indi2.F;
                                                    }
                                                    else
                                                    {
                                                        if (Elite.x[jj] + rho / 2 < 1)
                                                            indi2.x[jj] = Elite.x[jj] + rho / 2;
                                                        else
                                                            indi2.x[jj] = Elite.x[jj];
                                                        if (OF < 16 || OF == 22)
                                                            indi2.F = NFC.compute_fitness(indi2.x, OF, N);
                                                        else
                                                            indi2.F = NFC.compute_fitness(indi2.x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                        it++;
                                                        if (indi2.F >= Indi.F)
                                                        {
                                                            Indi.x = (double[])indi2.x.Clone();
                                                            Indi.F = indi2.F;
                                                        }
                                                    }
                                                }
                                                if (Indi.F >= Elite.F)
                                                {
                                                    Elite.x = (double[])Indi.x.Clone();
                                                    Elite.F = Indi.F;
                                                    updated = true;
                                                }
                                                else
                                                    rho /= 2;

                                            }
                                            #endregion
                                            if (it < NPop * Iteration)
                                            {
                                                if (updated)
                                                {
                                                    #region Middle_distance
                                                    for (int I = 0; I < k * N; I++)
                                                    {
                                                        Indi.x = NFC.Initialization_Hypercube(Elite.x, N, Rand, delta);
                                                        int i = Rand.Next(N);
                                                        Indi.x[i] = Elite.x[i];
                                                        int checker = 0;
                                                        while (Rand.NextDouble() <= Crprim)
                                                        {
                                                            Indi.x[i] = Elite.x[i];
                                                            i++;
                                                            if (i == N)
                                                                i = 0;
                                                            if (checker >= N)
                                                                break;
                                                            checker++;
                                                        }
                                                        //Indi.x = NFC.Node_treatment(Indi.x, N, Rand);
                                                        Solution_Treatment(boundary_method, ref Indi.x, Rand);
                                                        if (OF < 16 || OF == 22)
                                                            Indi.F = NFC.compute_fitness(Indi.x, OF, N);
                                                        else
                                                            Indi.F = NFC.compute_fitness(Indi.x, OF, N, Perm, Orthogonal_Matrix, M, o);

                                                        if (Indi.F > Elite.F)
                                                            updated = true;
                                                        if (Indi.F >= Elite.F)
                                                        {
                                                            Elite.x = (double[])Indi.x.Clone();
                                                            Elite.F = Indi.F;
                                                        }
                                                        it++;
                                                    }
                                                    #endregion
                                                }
                                                else
                                                {
                                                    #region Long_distance
                                                    Indi.x = NFC.node_constructing(N, Rand);
                                                    int i = Rand.Next(N);
                                                    Indi.x[i] = Elite.x[i];
                                                    int checker = 0;
                                                    while (Rand.NextDouble() <= Cr)
                                                    {
                                                        Indi.x[i] = Elite.x[i];
                                                        i++;
                                                        if (i == N)
                                                            i = 0;
                                                        if (checker >= N)
                                                            break;
                                                        checker++;
                                                    }
                                                    if (OF < 16 || OF == 22)
                                                        Indi.F = NFC.compute_fitness(Indi.x, OF, N);
                                                    else
                                                        Indi.F = NFC.compute_fitness(Indi.x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                    if (Indi.F >= Elite.F)
                                                    {
                                                        Elite.x = (double[])Indi.x.Clone();
                                                        Elite.F = Indi.F;
                                                        updated = true;
                                                    }
                                                    //////////////////////Function Evaluation increasing///////////////////
                                                    it++;
                                                    #endregion
                                                }
                                            }
                                        }
                                        costs[tt] = Elite.F;
                                    }
                                    sw = new StreamWriter(result);
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        sw.Write(costs[tt].ToString() + " ");
                                    }
                                    sw.Close();

                                }
                                //}
                                #endregion


                            }
                            #endregion
                            #region BBO
                            else if (method == "BBO")
                            {
                                double[] MutationRate = new double[NPop];
                                bool probflag = true;
                                double pmodify = 1;
                                int keep = 1;
                                //double[] dts = { 0.5, 2, 5, 50, 200 };
                                ////double mu_rate = 0.001;
                                //double[] Mus = { 0, 0.001, 0.005, 0.01, 0.1 };
                                double[] dts = { 200, 200, 200, 50, 5, 2 };
                                //double mu_rate = 0.001;
                                double[] Mus = { 0.005, 0.01, 0.001, 0.01, 0.01, 0.001 };
                                // double[]Dt={
                                double MaxParValue = 1; double MinParValue = -1;
                                double[,] chromkeep = new double[keep, N];
                                double[] fikeep = new double[keep];
                                int I = 1; int E = 1; int P = NPop;
                                //for(int ij=0;ij<5;ij++)
                                //    for (int iji = 0; iji < 5; iji++)
                                //    {
                                //double dt = dts[OF - 22];
                                //double mu_rate = Mus[OF - 22];
                                double dt = dts[0];
                                double mu_rate = Mus[0];
                                //double dt=Dt[
                                //ss2 = "Results\\fit_N_" + N.ToString() + "_bound_" + boundary_method.ToString() + "pop_" + NPop.ToString() + "_MU_" + mu_rate.ToString() + "_Dt" + dt.ToString() + method + "_test_" + Ntest.ToString() + "OF_" + OF.ToString() + ".txt";
                                #region file_existency
                                if (!File.Exists(result))
                                {
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        BBO[] BBO = new BBO[NPop];
                                        BBO[] Island = new OMOA_Namespace.BBO[NPop];
                                        BBO Best_indi = new BBO();
                                        Best_indi.F = Double.NegativeInfinity;
                                        Fmax = Double.NegativeInfinity; Fmin = 10000000000;
                                        for (int ii = 0; ii < NPop; ii++)
                                        {

                                            BBO[ii].x = NFC.node_constructing(N, Rand);
                                            if (OF < 16 || OF == 22)
                                                BBO[ii].F = NFC.compute_fitness(BBO[ii].x, OF, N);
                                            else
                                                BBO[ii].F = NFC.compute_fitness(BBO[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                            if (Fmax < BBO[ii].F)
                                                Fmax = BBO[ii].F;
                                        }
                                        BBO = NFC.mergesort(BBO, NPop);
                                        double[] Prob = new double[NPop];
                                        bool prob_flag = true;
                                        for (int ii = 0; ii < NPop; ii++)
                                        {
                                            Prob[ii] = 1 / BBO[ii].F;
                                        }
                                        for (int it = 0; it < Iteration; it++)
                                        {
                                            for (int ii = 0; ii < keep; ii++)
                                            {
                                                for (int jj = 0; jj < N; jj++)
                                                    chromkeep[ii, jj] = BBO[ii].x[jj];
                                                fikeep[ii] = BBO[ii].F;
                                            }
                                            NFC.GetSpeciesCounts(ref BBO);
                                            LS temp = NFC.GetLambdaMu(BBO, I, E);
                                            double[] lambda = (double[])temp.bias.Clone();
                                            double[] mu = (double[])temp.x.Clone();
                                            double sum = 0;
                                            if (prob_flag)
                                            {
                                                double[] probdot = new double[NPop];
                                                for (int jj = 0; jj < NPop; jj++)
                                                {
                                                    double ProbMinus = -1; double ProbPlus = -1;
                                                    double lambdaMinus = I * (1 - (BBO[jj].Species - 1) / NPop);
                                                    double muPlus = E * (BBO[jj].Species + 1) / NPop;
                                                    if (jj < NPop - 1)
                                                        ProbMinus = Prob[jj + 1];
                                                    else
                                                        ProbMinus = 0;
                                                    if (jj > 0)
                                                        ProbPlus = Prob[jj - 1];
                                                    else
                                                        ProbPlus = 0;
                                                    probdot[jj] = -(lambda[jj] + mu[jj]) * Prob[jj] + lambdaMinus * ProbMinus + muPlus * ProbPlus;
                                                    Prob[jj] = Prob[jj] + probdot[jj] * dt;
                                                    if (Prob[jj] < 0)
                                                        Prob[jj] = 0;
                                                    sum += Prob[jj];
                                                }
                                            }
                                            for (int ii = 0; ii < NPop; ii++)
                                                Prob[ii] = Prob[ii] / sum;

                                            double lambdaMin = 1000;
                                            double lambdaMax = -10000;
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                if (lambda[ii] < lambdaMin)
                                                    lambdaMin = lambda[ii];
                                                else if (lambda[ii] > lambdaMax)
                                                    lambdaMax = lambda[ii];
                                            }
                                            double lambdaScale = 0; double lambdaLower = 0.0; double lambdaUpper = 1;
                                            for (int k = 0; k < NPop; k++)
                                            {
                                                if (lambdaLower < lambda[k])
                                                    lambdaLower = lambda[k];
                                                else if (lambdaUpper > lambda[k])
                                                    lambdaUpper = lambda[k];
                                                if (Rand.NextDouble() > pmodify)
                                                    continue;
                                                lambdaScale = lambdaLower + (lambdaUpper - lambdaLower) * (lambda[k] - lambdaMin) / (lambdaMax - lambdaMin);

                                                //% Probabilistically input new information into habitat i
                                                Island[k].x = new double[N];
                                                for (int j = 0; j < N; j++)
                                                {
                                                    if (Rand.NextDouble() < lambdaScale)
                                                    {
                                                        //                % Pick a habitat from which to obtain a feature
                                                        sum = 0;
                                                        for (int ii = 0; ii < NPop; ii++)
                                                            sum += mu[ii];
                                                        double RandomNum = Rand.NextDouble() * sum;
                                                        double Select = mu[1];
                                                        int SelectIndex = 0;
                                                        while ((RandomNum > Select) && (SelectIndex < NPop - 1))
                                                        {
                                                            SelectIndex = SelectIndex + 1;
                                                            Select = Select + mu[SelectIndex];
                                                        }

                                                        Island[k].x[j] = BBO[SelectIndex].x[j];
                                                    }
                                                    else
                                                        Island[k].x[j] = BBO[k].x[j];
                                                }
                                            }

                                            if (probflag)
                                            {
                                                //        % Mutation
                                                double Pmax = -10000;
                                                for (int ii = 0; ii < NPop; ii++)
                                                {
                                                    if (Pmax < Prob[ii])
                                                        Pmax = Prob[ii];
                                                }
                                                for (int ii = 0; ii < NPop; ii++)
                                                {
                                                    MutationRate[ii] = mu_rate * (1 - Prob[ii] / Pmax);
                                                }
                                                //        % Mutate only the worst half of the solutions
                                                BBO = NFC.mergesort(BBO, NPop);
                                                int k = NPop / 2;
                                                for (; k < NPop; k++)
                                                    for (int parnum = 0; parnum < NPop; parnum++)
                                                        if (MutationRate[k] > Rand.NextDouble())
                                                            Island[k].x[parnum] = MinParValue + (MaxParValue - MinParValue + 1) * Rand.NextDouble();
                                            }
                                            //    % Replace the habitats with their new versions.
                                            for (int k = 0; k < NPop; k++)
                                            {
                                                for (int jj = 0; jj < N; jj++)
                                                {
                                                    BBO[k].x[jj] = Island[k].x[jj];
                                                }
                                                Solution_Treatment(boundary_method, ref BBO[k].x, Rand);
                                            }
                                            //    % Make sure each individual is legal.
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                if (OF < 16 || OF == 22)
                                                    BBO[ii].F = NFC.compute_fitness(BBO[ii].x, OF, N);
                                                else
                                                    BBO[ii].F = NFC.compute_fitness(BBO[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                if (Fmax < BBO[ii].F)
                                                    Fmax = BBO[ii].F;
                                            }
                                            //    % Calculate cost
                                            //    % Sort from best to worst
                                            BBO = NFC.mergesort(BBO, NPop);
                                            //    % Replace the worst with the previous generation's elites.

                                            for (int k = 0; k < keep; k++)
                                            {
                                                for (int jj = 0; jj < N; jj++)
                                                    BBO[NPop - (k + 1)].x[jj] = chromkeep[k, jj];
                                                BBO[NPop - (k + 1)].F = fikeep[k];
                                            }
                                            //    % Make sure the population does not have duplicates. 
                                            //    Population = ClearDups(Population, MaxParValue, MinParValue);
                                            //    % Sort from best to worst
                                            BBO = NFC.mergesort(BBO, NPop);
                                        }
                                        costs[tt] = BBO[0].F;
                                    }
                                    sw = new StreamWriter(result);
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        sw.Write(costs[tt].ToString() + " ");
                                    }
                                    sw.Close();

                                    //}
                                }
                                #endregion



                            }
                            #endregion
                            #region JADE
                            else if (method == "JADE")
                            {

                                double[] c_a = { 0.1, 0.1, 0.0001, 0, 0.0001, 0, 0, 0.1, 0, 0.0001, 0, 0.0001, 0.0001, 0, 0, 0.001, 0, 0.001, 0.0001, 0.01, 0, 0.001 };
                                double[] p_a = { 0.01, 0.3, 0.3, 0.3, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.3, 0.05, 0.05, 0, 0.05, 0.5, 0.05, 0.1, 0.5, 0.5, 0.05, 0.1 };//greediness

                                //double[] Fs = { 0.1, 0.5, 0.5, 0.1, 0.1, 0.1, 0.5, 0.5, 0.1, 0.5, 0.1, 0.1, 0.1,0, 0.1, 0.1, 0.1, 0.5, 2.0, 0.1, 0.1, 0.1 };
                                //double[] C_s = { 0.2, 0.2, 0.2, 0.4, 0.2, 0.4, 0.2, 0.2, 0.4, 0.2, 0.2, 0.2, 0.2, 0,0.4, 1, 1, 0.2, 0.2, 1, 0.4, 0.8 };
                                double[] Cr_A = new double[NPop];
                                double[] F_A = new double[NPop];
                                double Cr = 0.5;
                                double F = 0.5;
                                //double[] c_a = { 0, 0.0001, 0.001, 0.01, 0.1 };
                                //double[] p_a = { 0.01, 0.05, 0.1, 0.3, 0.5 };
                                //double[] c_a = { 0.1, 0, 0, 0.1, 0.001, 0.1 };
                                //double[] p_a = { 0.01, 0.05, 0.05, 0.05, 0.1, 0.05 };
                                //for (int CC = 0; CC < 5; CC++)
                                //    for (int PP = 0; PP < 5; PP++)
                                //    {
                                double c = c_a[OF ];
                                double p = p_a[OF ];
                                int p_group = Convert.ToInt32(100 * p % NPop);
                                //CD.Scale=0.1;
                                //if (c == 0.000001)
                                //    ss2 = "fit_N_" + N.ToString() + "_Pop_" + NPop.ToString() + "_c_" + 6.ToString() + "_p_" + p.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "OF_" + OF.ToString() + ".txt";
                                //else 
                                //ss2 = "Results\\fit_N_" + N.ToString() + "_bound_" + boundary_method.ToString() + "_Pop_" + NPop.ToString() + "_c_" + c.ToString() + "_p_" + p.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "OF_" + OF.ToString() + ".txt";
                                if (!File.Exists(result))
                                {
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        Cr_A = new double[NPop];
                                        F_A = new double[NPop];
                                        Cr = 0.5;
                                        F = 0.5;
                                        NF[] Nodes = new NF[NPop];
                                        NF[] V = new NF[NPop];
                                        NF[] U = new NF[NPop];
                                        NF[] Best = new NF[NPop];
                                        ArrayList SF = new ArrayList();
                                        ArrayList Scr = new ArrayList();
                                        ArrayList A_x = new ArrayList();
                                        ArrayList A_f = new ArrayList();
                                        NF Best_nodes = new NF();
                                        double[] distances = new double[NPop];
                                        Best_nodes.F = -1E250;
                                        for (int ii = 0; ii < NPop; ii++)
                                        {
                                            Nodes[ii].x = NFC.node_constructing(N, Rand);
                                            A_x.Add(Nodes[ii].x.Clone());
                                            if (OF < 16 || OF == 22)
                                                Nodes[ii].F = NFC.compute_fitness(Nodes[ii].x, OF, N);
                                            else
                                                Nodes[ii].F = NFC.compute_fitness(Nodes[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                            A_f.Add(Nodes[ii].F);
                                            /////////init V
                                            V[ii].x = new double[N];
                                            U[ii].x = new double[N];
                                        }
                                        ///////////////////////////////////////////////////////////////////////////////////////
                                        for (int it = 0; it < Iteration; it++)
                                        {
                                            SF.Clear();
                                            Scr.Clear();
                                            Best = NFC.mergesort(Nodes, NPop);
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                F_A[ii] = NFC.Cauchy(Rand, F, 0.1);
                                                Cr_A[ii] = NFC.Normal_Distribution(Rand, Cr, 0.1);
                                                int best_index = Rand.Next(p_group);
                                                NF X_best = Best[best_index];
                                                int index_r1 = Rand.Next(NPop);
                                                if (X_best.F == Nodes[index_r1].F)
                                                {
                                                    while (NFC.similarity(X_best.x, Nodes[index_r1].x))
                                                        index_r1 = Rand.Next(NPop);
                                                }
                                                int index_r2 = Rand.Next(A_f.Count);
                                                double Fit_A = Convert.ToDouble(A_f[index_r2]);
                                                if (X_best.F == Fit_A || Fit_A == Nodes[index_r1].F)
                                                {
                                                    double[] X = (double[])A_x[index_r2];
                                                    while (NFC.similarity(X_best.x, X))
                                                    {
                                                        index_r2 = Rand.Next(NPop);
                                                        X = (double[])A_x[index_r2];
                                                    }
                                                    while (NFC.similarity(Nodes[index_r1].x, X))
                                                    {
                                                        index_r2 = Rand.Next(NPop);
                                                        X = (double[])A_x[index_r2];
                                                    }
                                                }
                                                int J_rand = Rand.Next(N);
                                                for (int jj = 0; jj < N; jj++)
                                                {
                                                    double[] x = (double[])A_x[index_r2];
                                                    V[ii].x[jj] = Nodes[ii].x[jj] + F_A[ii] * (X_best.x[jj] - Nodes[ii].x[jj]) + F_A[ii] * (Nodes[index_r1].x[jj] - x[jj]);
                                                    if ((jj == J_rand) || (Rand.NextDouble() < Cr_A[ii]))
                                                        U[ii].x[jj] = V[ii].x[jj];
                                                    else
                                                        U[ii].x[jj] = Nodes[ii].x[jj];

                                                }
                                                //U[ii].x = NFC.Node_treatment(U[ii].x, N, Rand);
                                                Solution_Treatment(boundary_method, ref U[ii].x, Rand);
                                                if (OF < 16 || OF == 22)
                                                    U[ii].F = NFC.compute_fitness(U[ii].x, OF, N);
                                                else
                                                    U[ii].F = NFC.compute_fitness(U[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                if (U[ii].F >= Nodes[ii].F)
                                                {
                                                    int rand = Rand.Next(NPop);
                                                    A_x.Add((double[])Nodes[ii].x.Clone());
                                                    A_f.Add(Nodes[ii].F);
                                                    Nodes[ii].x = (double[])U[ii].x.Clone();
                                                    Nodes[ii].F = U[ii].F;
                                                    SF.Add(F_A[ii]);
                                                    Scr.Add(Cr_A[ii]);
                                                }
                                                if (Best_nodes.F < Nodes[ii].F)
                                                {
                                                    Best_nodes.F = Nodes[ii].F;
                                                    Best_nodes.x = (double[])Nodes[ii].x.Clone();

                                                }
                                            }

                                            int Overhead = A_x.Count - NPop;
                                            for (int ii = 0; ii < Overhead; ii++)
                                            {
                                                int rand = Rand.Next(A_x.Count);
                                                A_x.RemoveAt(rand);
                                                A_f.RemoveAt(rand);
                                            }

                                            Cr = (1 - c) * Cr + c * NFC.mean_simple(Scr);
                                            F = (1 - c) * F + c * NFC.mean_Lehmer(SF);
                                        }

                                        costs[tt] = Best_nodes.F;
                                        Range.Refresh();
                                    }

                                    sw = new StreamWriter(result);
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        sw.Write(costs[tt].ToString() + " ");
                                    }
                                    sw.Close();
                                }
                                //}
                            }

                            #endregion
                            #region PMS
                            if (method == "PMS")
                            {

                                double[] h = new double[N];
                                for (int ii = 0; ii < N; ii++)
                                    h[ii] = 0.1;
                                double[,] A = new double[N, N];
                                //double[] inheritance_factors = { 0.01, 0.05, 0.1, 0.3, 0.6 };
                                //double[] Alphas = {2,30,200 };
                                //double[] Bethas = { 0.7,1,5,10,50 };
                                bool updated = false;
                                //double[] rho_s = { 0.001, 0.05, 0.01, 0.05, 0.001, 0.01 };
                                //double[] inh_rate = { 0.3, 0.1, 0.1, 0.1, 0.1, 0.3 };
                                double[] rho_s = { 0.01, 0.01, 0.001, 0.001, 0.01, 0.1, 0.7, 0.01, 0.01, 0.01, 0.01, 0.01, 0.001, 0, 0.001, 0.01, 0.1, 0.001, 0.001, 0.01, 0.01, 0.1 };
                                double[] inh_rate = { 0.05, 0.01, 0.05, 0.001, 0.001, 0.05, 0.01, 0.001, 0.01, 0.05, 0.3, 0.05, 0.05, 0, 0.3, 0.01, 0.05, 0.001, 0.3, 0.05, 0.1, 0.01 };
                                //for (int rh = 0; rh < 5; rh++)
                                //    for (int cr = 0; cr < 5; cr++)
                                //    {
                                double inheritance_factor = inh_rate[OF ];
                                double root = N * inheritance_factor;
                                double Cr = 1.0 / (Math.Pow(2, 1.0 / root));
                                root = N * (1 - inheritance_factor);
                                int local_budget = 150;
                                double alpha = 2;
                                double betha = 0.5;
                                double rho = rho_s[OF ];
                                //ss2 = "Results\\fit_N_" + N.ToString() + "_bound_" + boundary_method.ToString() +"pop_" + NPop.ToString() + method + "inh_rate_" + inheritance_factor.ToString() + "rho_" + rho.ToString() + "_test_" + Ntest.ToString() + "OF_" + OF.ToString() + ".txt";
                                #region file_existency
                                if (!File.Exists(result))
                                {
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        rho = rho_s[OF];
                                        for (int ii = 0; ii < N; ii++)
                                            A[ii, ii] = 1;
                                        NF Elite = new NF();
                                        NF Indi = new NF();
                                        ////////////////////////////////////
                                        //////////////////////////////
                                        Elite.x = NFC.node_constructing(N, Rand);
                                        if (OF < 16 || OF == 22)
                                            Elite.F = NFC.compute_fitness(Elite.x, OF, N);
                                        else
                                            Elite.F = NFC.compute_fitness(Elite.x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                        for (int it = 0; it < Iteration * NPop; )
                                        {
                                            ///////////////////Long distance/////////////////////////////////
                                            updated = false;
                                            int global_budget_spent_gen = Convert.ToInt32(0.05 * Iteration * NPop);
                                            int gg = 0;
                                            while (!updated && it < Iteration * NPop && gg < global_budget_spent_gen)
                                            {
                                                #region Long_distance
                                                int checker = 0;
                                                Indi.x = NFC.node_constructing(N, Rand);
                                                int i = Rand.Next(N);
                                                Indi.x[i] = Elite.x[i];
                                                while (Rand.NextDouble() <= Cr)
                                                {
                                                    Indi.x[i] = Elite.x[i];
                                                    i++;
                                                    if (i == N)
                                                        i = 0;
                                                    if (checker >= N)
                                                        break;
                                                    checker++;
                                                }
                                                Solution_Treatment(boundary_method,ref Indi.x, Rand);
                                                if (OF < 16 || OF == 22)
                                                    Indi.F = NFC.compute_fitness(Indi.x, OF, N);
                                                else
                                                    Indi.F = NFC.compute_fitness(Indi.x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                if (Indi.F >= Elite.F)
                                                {
                                                    Elite.x = (double[])Indi.x.Clone();
                                                    Elite.F = Indi.F;
                                                    updated = true;
                                                }
                                                //////////////////////Function Evaluation increasing///////////////////
                                                it++;
                                                gg++;
                                                #endregion
                                            }
                                            /////////////////////////Short Distance////////////////////////////////////////
                                            double decisional_par = Rand.NextDouble();
                                            if (decisional_par < 0.5)
                                            {
                                                #region Short_region
                                                NF indi2 = new NF();
                                                for (int ii = 0; ii < local_budget; ii++)
                                                {
                                                    Indi.x = (double[])Elite.x.Clone();
                                                    indi2.x = (double[])Elite.x.Clone();
                                                    for (int jj = 0; jj < N; jj++)
                                                    {
                                                        if (Elite.x[jj] - rho < -1)
                                                            indi2.x[jj] = Elite.x[jj] - rho;
                                                        else
                                                            indi2.x[jj] = Elite.x[jj];
                                                        Solution_Treatment(boundary_method, ref indi2.x, Rand);
                                                        if (OF < 16 || OF == 22)
                                                            indi2.F = NFC.compute_fitness(indi2.x, OF, N);
                                                        else
                                                            indi2.F = NFC.compute_fitness(indi2.x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                        it++;
                                                        if (indi2.F >= Indi.F)
                                                        {
                                                            Indi.x = (double[])indi2.x.Clone();
                                                            Indi.F = indi2.F;
                                                        }
                                                        else
                                                        {
                                                            if (Elite.x[jj] + rho / 2 > 1)
                                                                indi2.x[jj] = Elite.x[jj] + rho / 2;
                                                            else
                                                                indi2.x[jj] = Elite.x[jj];
                                                            if (OF < 16 || OF == 22)
                                                                indi2.F = NFC.compute_fitness(indi2.x, OF, N);
                                                            else
                                                                indi2.F = NFC.compute_fitness(indi2.x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                            it++;
                                                            if (indi2.F >= Indi.F)
                                                            {
                                                                Indi.x = (double[])indi2.x.Clone();
                                                                Indi.F = indi2.F;
                                                            }
                                                        }
                                                    }
                                                    if (Indi.F >= Elite.F)
                                                    {
                                                        Elite.x = (double[])Indi.x.Clone();
                                                        Elite.F = Indi.F;
                                                    }
                                                    else
                                                        rho /= 2;

                                                }
                                                #endregion
                                            }
                                            else
                                            {
                                                NFC.rosenbrock(50000, alpha, betha, N, ref Elite.F, ref Elite.x, Rand, OF, Perm, Orthogonal_Matrix, M, o, ref it);
                                            }
                                        }
                                        costs[tt] = Elite.F;
                                    }
                                    sw = new StreamWriter(result);
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        sw.Write(costs[tt].ToString() + " ");
                                    }
                                    sw.Close();

                                    //}
                                }
                                #endregion


                            }
                            #endregion
                            #region CCPSO2
                            else if (method == "CCPSO2")
                            {
                                //StreamReader SR = new StreamReader("cauchy.txt");
                                //string all = SR.ReadToEnd();
                                //string[] all_with_space = all.Split('\n', ',');
                                //double[] Cauchy = new double[50000];
                                //for (int ii = 0; ii < 50000; ii++)
                                //{
                                //    Cauchy[ii] = Convert.ToDouble(all_with_space[ii]);
                                //}
                                double[,] S1_A = { { 0.02, 0.04, 0.08, 0.16, 0.32 }, { 0.02, 0.05, 0.1, 0.5, 1 }, { 0.04, 0.16, 0.32, 0.64, 1 }, { 0.1, 0.3, 0.5, 0.7, 0.9 }, { 0.02, 0.2, 0.4, 0.8, 1 } };
                                int[] S = { 2, 5, 10, 50, 100 };
                                int[] SSS = { 2, 4, 4, 1, 2, 4 };
                                //double[] p_s = { 0, 0.001, 0.1, 0.5, 1 };
                                double[] p_s = { 1, 0.5, 1, 1, 1, 1 };
                                //for (int sss = 0; sss < 5; sss++)
                                //{
                                for (int ii = 0; ii < S.Length; ii++)
                                    S[ii] = Convert.ToInt32(S1_A[SSS[1], ii] * N);
                                //for (int pp = 0; pp < 5; pp++)
                                //{
                                double p = p_s[0];// it has a chance to be a choice as a parameter
                                NF[] Pop = new NF[NPop];
                                int[] Indexes = new int[N];
                                //ss2 = "Results\\fit_N_" + N.ToString() + "_bound_" + boundary_method.ToString() + "_Pop_" + NPop.ToString() + "p_" + p.ToString() + "S" + SSS[1].ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "OF_" + OF.ToString() + ".txt";
                                if (!File.Exists(result))
                                {
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        int ind = Rand.Next(S.Length);
                                        int s = S[ind];
                                        int K = N / s;
                                        Swarms[] SWs = new Swarms[K];
                                        Swarms[] LocalBest = new Swarms[K];
                                        Swarms[] PersonalBest = new Swarms[K];
                                        NF[] PersonBest = new NF[NPop];
                                        NF Ybest = new NF();
                                        Ybest.F = -1E290;
                                        Swarms[] PYbest = new Swarms[K];
                                        for (int ii = 0; ii < NPop; ii++)
                                        {
                                            Pop[ii].x = NFC.node_constructing(N, Rand);
                                            PersonBest[ii].x = (double[])Pop[ii].x.Clone();
                                            if (OF < 16 || OF == 22)
                                                Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N);
                                            else
                                                Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                            if (Ybest.F < Pop[ii].F)
                                            {
                                                Ybest.x = (double[])Pop[ii].x.Clone();
                                                Ybest.F = Pop[ii].F;
                                            }
                                            PersonBest[ii].F = Pop[ii].F;
                                        }
                                        double p_yfirst = Ybest.F;
                                        for (int it = 0; it < NPop * Iteration; )
                                        {
                                            if (p_yfirst >= Ybest.F && it != 0)
                                            {
                                                ind = Rand.Next(S.Length);
                                                s = S[ind];
                                                K = N / s;
                                                SWs = new Swarms[K];
                                                LocalBest = new Swarms[K];
                                                PersonalBest = new Swarms[K];
                                                PYbest = new Swarms[K];
                                            }
                                            p_yfirst = Ybest.F;
                                            Indexes = NFC.sq_rand_gen2(N, N, Rand);
                                            int ll = 0;
                                            for (int kk = 0; kk < K; kk++)
                                            {
                                                SWs[kk].particles = new NF[NPop];
                                                LocalBest[kk].particles = new NF[NPop];
                                                PersonalBest[kk].particles = new NF[NPop];
                                                PYbest[kk].particles = new NF[1];
                                                for (int jj = 0; jj < NPop; jj++)
                                                {
                                                    SWs[kk].particles[jj].x = new double[s];
                                                    LocalBest[kk].particles[jj].x = new double[s];
                                                    PersonalBest[kk].particles[jj].x = new double[s];
                                                    PYbest[kk].particles[0].x = new double[s];
                                                    for (int ii = 0; ii < s; ii++)
                                                    {
                                                        SWs[kk].particles[jj].x[ii] = Pop[jj].x[Indexes[kk * s + ii]];
                                                        PersonalBest[kk].particles[jj].x[ii] = PersonBest[jj].x[Indexes[kk * s + ii]];
                                                    }
                                                    Solution_Treatment(boundary_method, ref PersonalBest[kk].particles[jj].x, Rand);
                                                }
                                                Solution_Treatment(boundary_method, ref Ybest.x, Rand);
                                                for (int ii = 0; ii < s; ii++)
                                                {
                                                    PYbest[kk].particles[0].x[ii] = Ybest.x[ll++];
                                                }

                                                PYbest[kk].particles[0].F = Ybest.F;
                                            }
                                            for (int ss = 0; ss < K; ss++)
                                            {
                                                for (int jj = 0; jj < NPop; jj++)
                                                {
                                                    double[] constructing = NFC.Concatenation(Ybest.x, N, ss * s, s, SWs[ss].particles[jj].x);
                                                    double[] personal_best = NFC.Concatenation(Ybest.x, N, ss * s, s, PersonalBest[ss].particles[jj].x);
                                                    double[] y_best = NFC.Concatenation(Ybest.x, N, ss * s, s, PYbest[ss].particles[0].x);
                                                    if (OF < 16 || OF == 22)
                                                        SWs[ss].particles[jj].F = NFC.compute_fitness(constructing, OF, N);
                                                    else
                                                        SWs[ss].particles[jj].F = NFC.compute_fitness(constructing, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                    if (OF < 16 || OF == 22)
                                                        PersonalBest[ss].particles[jj].F = NFC.compute_fitness(personal_best, OF, N);
                                                    else
                                                        PersonalBest[ss].particles[jj].F = NFC.compute_fitness(personal_best, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                    if (PersonalBest[ss].particles[jj].F < SWs[ss].particles[jj].F)
                                                    {
                                                        PersonalBest[ss].particles[jj].F = SWs[ss].particles[jj].F;
                                                        PersonalBest[ss].particles[jj].x = (double[])SWs[ss].particles[jj].x.Clone();
                                                    }
                                                    if (OF < 16 || OF == 22)
                                                        PYbest[ss].particles[0].F = NFC.compute_fitness(y_best, OF, N);
                                                    else
                                                        PYbest[ss].particles[0].F = NFC.compute_fitness(y_best, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                    if (PersonalBest[ss].particles[jj].F > PYbest[ss].particles[0].F)
                                                    {
                                                        PYbest[ss].particles[0].F = PersonalBest[ss].particles[jj].F;
                                                        PYbest[ss].particles[0].x = (double[])PersonalBest[ss].particles[jj].x.Clone();
                                                    }
                                                    it = it + 3;
                                                }
                                                /////////////////////////////////////////////////
                                                for (int jj = 0; jj < NPop; jj++)
                                                {
                                                    int[] Nei = NFC.Ring_withmyself(jj, NPop);
                                                    double max = -1E250;
                                                    for (int ii = 0; ii < Nei.Length; ii++)
                                                    {
                                                        if (PersonalBest[ss].particles[Nei[ii]].F > max)
                                                        {
                                                            LocalBest[ss].particles[jj].F = PersonalBest[ss].particles[Nei[ii]].F;
                                                            LocalBest[ss].particles[jj].x = (double[])PersonalBest[ss].particles[Nei[ii]].x.Clone();
                                                            max = PersonalBest[ss].particles[Nei[ii]].F;
                                                        }
                                                    }
                                                }
                                                if (PYbest[ss].particles[0].F > Ybest.F)
                                                {
                                                    Ybest.F = PYbest[ss].particles[0].F;
                                                    int mm = 0;
                                                    for (int ii = ss * s; ii < ss * s + s; ii++)
                                                        Ybest.x[ii] = PYbest[ss].particles[0].x[mm++];
                                                }


                                            }
                                            for (int ss = 0; ss < K; ss++)
                                            {
                                                for (int jj = 0; jj < NPop; jj++)
                                                {
                                                    for (int kk = 0; kk < s; kk++)
                                                    {
                                                        if (Rand.NextDouble() <= p)
                                                            SWs[ss].particles[jj].x[kk] = PersonalBest[ss].particles[jj].x[kk] + NFC.Cauchy_mu(Rand) * (PersonalBest[ss].particles[jj].x[kk] - LocalBest[ss].particles[jj].x[kk]);
                                                        else
                                                            SWs[ss].particles[jj].x[kk] = LocalBest[ss].particles[jj].x[kk] + NFC.Normal_Distribution(Rand, 0, 1) * (PersonalBest[ss].particles[jj].x[kk] - LocalBest[ss].particles[jj].x[kk]);
                                                        Pop[jj].x[ss * s + kk] = SWs[ss].particles[jj].x[kk];
                                                        PersonBest[jj].x[ss * s + kk] = PersonalBest[ss].particles[jj].x[kk];
                                                    }
                                                }
                                            }
                                            for (int jj = 0; jj < NPop; jj++)
                                            {
                                                NFC.Node_treatment(Pop[jj].x, N, Rand);
                                            }

                                        }

                                        costs[tt] = Ybest.F;
                                    }
                                    sw = new StreamWriter(result);
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        sw.Write(costs[tt].ToString() + " ");
                                    }
                                    sw.Close();
                                }
                                //    }
                                //}
                            }

                            #endregion
                            #region MASSW
                            else if (method == "MASSW")
                            {
                                O = 50000;
                                //double[] N_LS_GL = { 0.7, 0.1, 0.5, 0.7, 0.7, 0.9 };
                                //double[] Mu_rate_s = { 0.05, 0.0, 0.0, 0.0, 0.01, 0.1 };
                                double[] N_LS_GL = { 0.3, 0.3, 0.5, 0.3, 0.3, 0.1, 0.3, 0.5, 0.3, 0.5, 0.1, 0.1, 0.3, 0, 0.3, 0.1, 0.5, 0.5, 0.5, 0.7, 0.3, 0.5 };
                                double[] Mu_rate_s = { 0.001, 0.1, 0.001, 0, 0, 0, 0.5, 0.001, 0.01, 0.01, 0.01, 0, 0.1, 0, 0.001, 0.01, 0.01, 0.001, 0.5, 0.1, 0.01, 0.001 };
                                //string[] Structs = { "Ring", "Star", "Btree", "Bipartie", "Cluster", "Grid", "Ladder", "Cross-Ladder", "Cellular" };
                                //for (int Portion_index = 0; Portion_index < 5; Portion_index++)
                                //    for (int M_index = 0; M_index < 5; M_index++)
                                //    {
                                //for (int SS = 0; SS < 9; SS++)
                                //{
                                //string Struct = Structs[SS];
                                //double[] N_LS_GL = { 0.5, 0.5, 0.9, 0.3, 0.5, 0.5, 0.1, 0.5, 0.3, 0.9, 0.9, 0.3, 0.5, 0, 0.3, 0.5, 0.5, 0.5, 0.1, 0.9, 0.3, 0.5 };
                                //double[] Mu_rate_s = { .15, 0.005, 0.005, 0.05, 0.05, 0.05, 0.05, 0.005, 0.005, 0.05, 0.05, 0.05, 0.01, 0, 0.05, 0.005, 0.01, 0.005, 0.05, 0.01, 0.05, 0.05 };
                                double LS_LG = N_LS_GL[OF ];
                                double Mu_rate = Mu_rate_s[OF ];
                                Iteration = Convert.ToInt32((O * (1 - LS_LG)) / NPop) + 1;
                                int I_str = Convert.ToInt32(O * LS_LG / Iteration);
                                int best_index = 0;
                                int local_search = 0;
                                // ss2 = "fit_N_" + N.ToString() + "_Mr_" + Mu_rate.ToString() + "LSLG" + LS_LG.ToString() + "_Pop_" + NPop.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString()+"Struct_"+Struct + "_OF" + OF.ToString() + ".txt";
                                //ss2 = "Results\\fit_N_" + N.ToString() + "_bound_" + boundary_method.ToString() + "_Mr_" + Mu_rate.ToString() + "LSLG" + LS_LG.ToString() + "_Pop_" + NPop.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_OF" + OF.ToString() + ".txt";
                                if (!File.Exists(result))
                                {
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        LS[] Pop = new LS[NPop];
                                        NF Best_nodes = new NF();
                                        NF[] Breed = new NF[NPop];
                                        //NF[] Breed = new NF[2*NPop];//for laplacian operator
                                        Best_nodes.F = -1E250;
                                        //LS[] LS_S = new LS[N];
                                        for (int ii = 0; ii < NPop; ii++)
                                        {
                                            Pop[ii].x = NFC.node_constructing(N, Rand);

                                            if (OF < 16 || OF == 22)
                                                Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N);
                                            else
                                                Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                            Pop[ii].local_search = true;

                                        }
                                        Pop = NFC.mergesort(Pop, NPop);
                                        for (int it = 0; it < Iteration; it++)
                                        {
                                            //perform Steady State GA
                                            #region generate Breeds
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                if (Mu_rate > Rand.NextDouble())
                                                    NFC.Mutation_BGA(ref Pop[ii], Pop, N, Rand);
                                                Solution_Treatment(boundary_method,ref Pop[ii].x, Rand);
                                                if (OF < 16 || OF == 22)
                                                    Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N);
                                                else
                                                    Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                            }
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                int parent1 = Rand.Next(NPop);
                                                Breed[ii].x = new double[N];
                                                //Breed[ii+1].x = new double[N];
                                                //int[] S = NFC.Structures(ii, NPop, Struct);
                                                int parent2 = NFC.negative_assortative_mating_strategy(Pop[parent1].x, Pop, Rand);
                                                Breed[ii].x = NFC.BLX_Crossover(Pop[parent1].x, Pop[parent2].x);
                                                //laplac temp = NFC.Laplacian_recombination(Pop[parent1], Pop[parent2], N, Rand, 0, 0.5);
                                                //Breed[ii].x = (double[])temp.par1.Clone();
                                                //Breed[ii+1].x = (double[])temp.par2.Clone();
                                                Breed[ii].x = NFC.Node_treatment(Breed[ii].x, N, Rand);
                                                //Breed[ii + 1].x = NFC.Node_treatment(Breed[ii + 1].x, N, Rand);
                                                if (OF < 16 || OF == 22)
                                                    Breed[ii].F = NFC.compute_fitness(Breed[ii].x, OF, N);
                                                else
                                                    Breed[ii].F = NFC.compute_fitness(Breed[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                //if (OF < 16 || OF == 22)
                                                //    Breed[ii+1].F = NFC.compute_fitness(Breed[ii+1].x, OF, N);
                                                //else
                                                //    Breed[ii+1].F = NFC.compute_fitness(Breed[ii+1].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                            }
                                            #endregion
                                            // Replacement Strategy

                                            Pop = NFC.mergesort(Pop, NPop);
                                            Breed = NFC.mergesort(Breed, NPop);
                                            if (Breed[0].F > Pop[NPop - 1].F)
                                            {
                                                Pop[NPop - 1].x = (double[])Breed[0].x.Clone();
                                                Pop[NPop - 1].F = Breed[0].F;
                                                Pop[NPop - 1].local_search = true;
                                            }
                                            //
                                            //Pick the best individual of Population
                                            if (local_search == NPop - 1)
                                            {
                                                for (int ii = 1; ii < NPop; ii++)
                                                {
                                                    Pop[ii].x = NFC.node_constructing(N, Rand);
                                                    if (OF < 16 || OF == 22)
                                                        Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N);
                                                    else
                                                        Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                    Pop[ii].local_search = true;
                                                    Pop[ii].visited = 0;
                                                }
                                                local_search = 0;
                                            }

                                            LS Best_Indi_BLS = new LS();
                                            Pop = NFC.mergesort(Pop, NPop);
                                            for (int ii = 0; ii < NPop; ii++)
                                                if (Pop[ii].local_search)
                                                {
                                                    Best_Indi_BLS.F = Pop[ii].F;
                                                    Best_Indi_BLS.x = (double[])Pop[ii].x.Clone();
                                                    Best_Indi_BLS.visited = Pop[ii].visited;
                                                    Best_Indi_BLS.rho = Pop[ii].rho;
                                                    if (Pop[ii].visited == 1)
                                                        Best_Indi_BLS.bias = (double[])Pop[ii].bias.Clone();
                                                    best_index = ii;
                                                    break;
                                                }
                                            // Initialize Bias and Rho
                                            if (Best_Indi_BLS.visited == 0)
                                            {
                                                Best_Indi_BLS.bias = new double[N];
                                                Best_Indi_BLS.rho = 0.01;
                                            }
                                            // Local Search Operator
                                            LS Best_Indi_ALS = new LS();
                                            Best_Indi_ALS.x = (double[])Best_Indi_BLS.x.Clone();
                                            Best_Indi_ALS.F = Best_Indi_BLS.F;
                                            Best_Indi_ALS.rho = Best_Indi_BLS.rho;
                                            Best_Indi_ALS.visited = Best_Indi_BLS.visited;
                                            Best_Indi_ALS.bias = (double[])Best_Indi_BLS.bias.Clone();
                                            NFC.Subgrouping_SolisWets(ref Best_Indi_ALS.x, ref Best_Indi_ALS.F, N, ref Best_Indi_ALS.bias, ref Best_Indi_ALS.rho, I_str, Rand, OF, Perm, Orthogonal_Matrix, M, o);

                                            // Remove Unchanged particle from population
                                            if (Best_Indi_BLS.F - Best_Indi_ALS.F == 0)
                                            {
                                                Pop[best_index].local_search = false;
                                                //  Pop[best_index].visited = 1;
                                                local_search++;
                                            }
                                            else
                                            {
                                                Pop[best_index].x = (double[])Best_Indi_ALS.x.Clone();
                                                Pop[best_index].visited = 1;
                                                Pop[best_index].F = Best_Indi_ALS.F;
                                                Pop[best_index].rho = Best_Indi_ALS.rho;
                                                Pop[best_index].bias = (double[])Best_Indi_ALS.bias.Clone();
                                                Pop[best_index].local_search = true;
                                            }

                                            Pop = NFC.mergesort(Pop, NPop);
                                            if (Best_nodes.F < Pop[0].F)
                                            {
                                                Best_nodes.F = Pop[0].F;
                                                Best_nodes.x = (double[])Pop[0].x.Clone();
                                            }

                                            /////////////////////

                                        }
                                        costs[tt] = Best_nodes.F;

                                    }


                                    sw = new StreamWriter(result);
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        sw.Write(costs[tt].ToString() + " ");
                                    }
                                    sw.Close();
                                }
                                //}

                            }

                            #endregion
                            #region MASSW_structure
                            else if (method == "MASSW_S")
                            {
                                O = NPop * Convert.ToInt32(NI.Text);
                                //double[] N_LS_GL = { 0.1, 0.3, 0.5, 0.7, 0.9 };
                                //double[] Mu_rate_s = {0.00001,0.001,0.01,0.1,0.5 };
                                double[] N_LS_GL = { 0.5, 0.7 };
                                double[] Mu_rate_s = { 0.01, 0.1, 0.5 };
                                string[] Structs = { "Ring", "Star", "Btree", "Bipartie", "Cluster", "Grid", "Ladder", "Cross-Ladder", "Cellular", "basic" };
                                //string[] Structs = { "Cellular" };
                                //for (int Portion_index = 0; Portion_index < 5; Portion_index++)
                                //    for (int M_index = 0; M_index < 5; M_index++)
                                //    {
                                for (int SS = 0; SS < Structs.Length; SS++)
                                {
                                    string Struct = Structs[SS];
                                    //double[] N_LS_GL = { 0.3,0.5,0.5,0.3,0.5,0.5,0.5,0.5,0.3,0.5,0.7,0.3,0.5,0,0.3,0.5,0.5,0.5,0.5,0.7,0.3,0.5 };
                                    //double[] Mu_rate_s = { 5, 0.5, 0.001, 0.001, 5, 5, 0.5, 0.5, 0.001, 0.01, 0.001, 0.01, 0.5, 0, 0.001, 0.1, 0.1, 0.5, 0.5, 5, 5, 0.1 };
                                    double LS_LG = 0.7;
                                    double Mu_rate = 0.5;
                                    Iteration = Convert.ToInt32((O * (1 - LS_LG)) / NPop) + 1;
                                    int I_str = Convert.ToInt32(O * LS_LG / Iteration);
                                    int best_index = 0;
                                    int local_search = 0;
                                    double[,] Diversity = new double[Ntest, Iteration];
                                    double[,] Fit = new double[Ntest, Iteration];
                                    //ss2 = "fit_N_" + N.ToString() + "_Mr_" + Mu_rate.ToString() + "LSLG" + LS_LG.ToString() + "_Pop_" + NPop.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "Struct_" + Struct + "_OF" + OF.ToString() + ".txt";
                                    //string ss = "Diver_N_" + N.ToString() + "_Mr_" + Mu_rate.ToString() + "LSLG" + LS_LG.ToString() + "_Pop_" + NPop.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "Struct_" + Struct + "_OF" + OF.ToString() + ".txt";
                                    //string ss1 = "Fitness_N_" + N.ToString() + "_Mr_" + Mu_rate.ToString() + "LSLG" + LS_LG.ToString() + "_Pop_" + NPop.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "Struct_" + Struct + "_OF" + OF.ToString() + ".txt";
                                    if (!File.Exists(result))
                                    {
                                        for (int tt = 0; tt < Ntest; tt++)
                                        {
                                            LS[] Pop = new LS[NPop];
                                            NF Best_nodes = new NF();
                                            NF Breed = new NF();
                                            Best_nodes.F = -1E250;
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                Pop[ii].x = NFC.node_constructing(N, Rand);

                                                if (OF < 16 || OF == 22)
                                                    Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N);
                                                else
                                                    Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                Pop[ii].local_search = true;

                                            }
                                            Pop = NFC.mergesort(Pop, NPop);
                                            for (int it = 0; it < Iteration; it++)
                                            {
                                                //perform Steady State GA
                                                //for (int ii = 0; ii < NPop; ii++)
                                                //{
                                                int parent1 = Rand.Next(NPop);
                                                int parent2 = -1;
                                                if (Struct != "basic")
                                                {
                                                    int[] S = NFC.Structures(parent1, NPop, Struct);
                                                    parent2 = NFC.negative_assortative_mating_strategy2(Pop[parent1].x, Pop, Rand, S);
                                                }
                                                else
                                                {
                                                    parent2 = NFC.negative_assortative_mating_strategy(Pop[parent1].x, Pop, Rand);
                                                }
                                                Breed.x = NFC.BLX_Crossover(Pop[parent1].x, Pop[parent2].x);
                                                if (Rand.NextDouble() < Mu_rate)
                                                    NFC.Mutation_BGA2(ref Breed, Pop, N, Rand);
                                                if (OF < 16 || OF == 22)
                                                    Breed.F = NFC.compute_fitness(Breed.x, OF, N);
                                                else
                                                    Breed.F = NFC.compute_fitness(Breed.x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                Breed.x = NFC.Node_treatment(Breed.x, N, Rand);
                                                Pop = NFC.mergesort(Pop, NPop);
                                                if (Breed.F > Pop[NPop - 1].F)
                                                {
                                                    Pop[NPop - 1].x = (double[])Breed.x.Clone();
                                                    Pop[NPop - 1].F = Breed.F;
                                                    Pop[NPop - 1].local_search = true;
                                                }
                                                //
                                                //Pick the best individual of Population
                                                if (local_search == NPop - 1)
                                                {
                                                    for (int ii = 1; ii < NPop; ii++)
                                                    {
                                                        Pop[ii].x = NFC.node_constructing(N, Rand);
                                                        if (OF < 16 || OF == 22)
                                                            Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N);
                                                        else
                                                            Pop[ii].F = NFC.compute_fitness(Pop[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                        Pop[ii].local_search = true;
                                                        Pop[ii].visited = 0;
                                                    }
                                                    local_search = 0;
                                                }

                                                LS Best_Indi_BLS = new LS();
                                                Pop = NFC.mergesort(Pop, NPop);
                                                for (int ii = 0; ii < NPop; ii++)
                                                    if (Pop[ii].local_search)
                                                    {
                                                        Best_Indi_BLS.F = Pop[ii].F;
                                                        Best_Indi_BLS.x = (double[])Pop[ii].x.Clone();
                                                        Best_Indi_BLS.visited = Pop[ii].visited;
                                                        Best_Indi_BLS.rho = Pop[ii].rho;
                                                        if (Pop[ii].visited == 1)
                                                            Best_Indi_BLS.bias = (double[])Pop[ii].bias.Clone();
                                                        best_index = ii;
                                                        break;
                                                    }
                                                // Initialize Bias and Rho
                                                if (Best_Indi_BLS.visited == 0)
                                                {
                                                    Best_Indi_BLS.bias = new double[N];
                                                    Best_Indi_BLS.rho = 0.01;
                                                }
                                                // Local Search Operator
                                                LS Best_Indi_ALS = new LS();
                                                Best_Indi_ALS.x = (double[])Best_Indi_BLS.x.Clone();
                                                Best_Indi_ALS.F = Best_Indi_BLS.F;
                                                Best_Indi_ALS.rho = Best_Indi_BLS.rho;
                                                Best_Indi_ALS.visited = Best_Indi_BLS.visited;
                                                Best_Indi_ALS.bias = (double[])Best_Indi_BLS.bias.Clone();
                                                NFC.Subgrouping_SolisWets(ref Best_Indi_ALS.x, ref Best_Indi_ALS.F, N, ref Best_Indi_ALS.bias, ref Best_Indi_ALS.rho, I_str, Rand, OF, Perm, Orthogonal_Matrix, M, o);

                                                // Remove Unchanged particle from population
                                                if (Best_Indi_BLS.F - Best_Indi_ALS.F == 0)
                                                {
                                                    Pop[best_index].local_search = false;
                                                    //  Pop[best_index].visited = 1;
                                                    local_search++;
                                                }
                                                else
                                                {
                                                    Pop[best_index].x = (double[])Best_Indi_ALS.x.Clone();
                                                    Pop[best_index].visited = 1;
                                                    Pop[best_index].F = Best_Indi_ALS.F;
                                                    Pop[best_index].rho = Best_Indi_ALS.rho;
                                                    Pop[best_index].bias = (double[])Best_Indi_ALS.bias.Clone();
                                                    Pop[best_index].local_search = true;
                                                }

                                                Pop = NFC.mergesort(Pop, NPop);
                                                if (Best_nodes.F < Pop[0].F)
                                                {
                                                    Best_nodes.F = Pop[0].F;
                                                    Best_nodes.x = (double[])Pop[0].x.Clone();
                                                }

                                                /////////////////////
                                                Diversity[tt, it] = NFC.diversity(Pop, NPop, N);
                                                Fit[tt, it] = Best_nodes.F;
                                            }
                                            costs[tt] = Best_nodes.F;
                                        }


                                        sw = new StreamWriter(result);
                                        for (int tt = 0; tt < Ntest; tt++)
                                        {
                                            sw.Write(costs[tt].ToString() + " ");
                                        }
                                        sw.Close();
                                        //sw = new StreamWriter(ss);
                                        //for (int tt = 0; tt < Ntest; tt++)
                                        //{
                                        //    for (int ii = 0; ii < Iteration; ii++)
                                        //    {
                                        //        sw.Write(Diversity[tt, ii].ToString() + " ");
                                        //    }
                                        //    sw.WriteLine();
                                        //}

                                        //sw.Close();
                                        sw = new StreamWriter(result);
                                        for (int tt = 0; tt < Ntest; tt++)
                                        {
                                            for (int ii = 0; ii < Iteration; ii++)
                                            {
                                                sw.Write(Fit[tt, ii].ToString() + " ");
                                            }
                                            sw.WriteLine();
                                        }
                                        sw.Close();
                                    }
                                    //}
                                }

                            }

                            #endregion
                            #region FMOA2013
                            else if (method == "FMOA2013")
                            {
                                //function_evaluation2013 class_new = new function_evaluation2013(ub, lb, dimension);

                                ////////////////////////////////////Initial parameter values/////////////////
                                //double Jr = Convert.ToDouble(Jumping_rate.Text);
                                //int acceleration = Convert.ToInt32(Ac.Text);
                                //int Distance = Convert.ToInt32(D.Text);
                                //double al1 = Convert.ToDouble(Al.Text);
                                //double al2 = Convert.ToDouble(Al2.Text);
                                //double al_step = Convert.ToDouble(Al_step.Text);
                                //double in1 = Convert.ToDouble(Inten.Text);
                                //double in2 = Convert.ToDouble(Inten2.Text);
                                //double in_step = Convert.ToDouble(Inten_step.Text);
                                //double short_range = 0.01;
                                //double[] Intinsities = { 0.001,0.01,0.2,0.6,1 };
                                //double[] Alphas = {0.001, 0.01,0.2,0.6,1 };
                                //double[] Jrs = { 0.3,0.45,0.6 };
                                //int ITT = 0;
                                //double[] Intinsities = { 0.6, 1.0, 1.0, 0.2, 0.6, 1.0, 0.2, 1.0, 1.0, 1.0, 1.0, 0.6, 1.0, 0, 1.0, 0.2, 0.2, 0.0010, 0.0010, 0.00100, 0.0010, 1.0 };
                                //double[] Alphas = { 0.0200000000000000, 0.400000000000000, 0.400000000000000, 0.200000000000000, 1, 0.800000000000000, 0.200000000000000, 0.200000000000000, 0.6, 0.4, 0.4, 1, 0.2, 0, 0.6, 1, 0.6, 0.4, 0.0010, 0.01000, 0.2000, 1 };
                                //double[] Jrs = { 0.6, 0.3, 0.45, 0.3, 0.3, 0.3, 0.45, 0.3, 0.3, 0.3, 0.3, 0.3, 0.45, 0, 0.3000, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.3 };
                                //int[] POP = { 5, 10, 25, 50, 100, 200, 250, 500 };
                                //int[] ITs = { 10000, 5000, 2000, 1000, 500, 250, 200, 100 };
                                //for (int al_index = 0; al_index < Alphas.Length; al_index++)
                                //    for (int int_index = 0; int_index < Intinsities.Length; int_index++)
                                //        for (int II = 0; II < 3; II++)
                                //        {
                                //double[] fitnesses = new double[Iteration];
                                //double intensity = 0.5;
                                //double alpha = 0.5;
                                //Jr = 0.3;
                                //double c = 0.5;
                                //double[] IntenI = new double[NPop];
                                //double[] AlphaI = new double[NPop];
                                //double[] JrI = new double[NPop];
                                //ArrayList IntenAr = new ArrayList();
                                //ArrayList AlphaAR = new ArrayList();
                                //ArrayList JrAR = new ArrayList();
                                //double[] PervFitness = new double[NPop];
                                //double[,] Js = new double[Ntest, Iteration];
                                //double[,] Rhos = new double[Ntest, Iteration];
                                //double[,] Alphas = new double[Ntest, Iteration];
                                //double[] ranges = new double[] { 0.001, 0.005,0.01,0.05 ,0.1 };
                                ////NPop = POP[II];
                                ////Iteration = ITs[II];
                                ////ss2 = "fit_N_" + N.ToString() + "_Pop_" + NPop.ToString() + "_intensity" + intensity.ToString() + "alpha_" + alpha.ToString() + "_c_" + c.ToString() + method + "_Acc_" + Ac.Text + "JR" + Jr.ToString() + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_struct_" + structure + "OF_" + OF.ToString() + ".txt";
                                //string ss1 = "JR_N_" + N.ToString() + "_Pop_" + NPop.ToString() + "_intensity" + intensity.ToString() + "alpha_" + alpha.ToString() + "_c_" + c.ToString() + method + "_Acc_" + Ac.Text + "JR" + Jr.ToString() + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_struct_" + structure + "OF_" + OF.ToString() + ".txt";
                                //string ss3 = "Rho_N_" + N.ToString() + "_Pop_" + NPop.ToString() + "_intensity" + intensity.ToString() + "alpha_" + alpha.ToString() + "_c_" + c.ToString() + method + "_Acc_" + Ac.Text + "JR" + Jr.ToString() + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_struct_" + structure + "OF_" + OF.ToString() + ".txt";
                                //string ss4 = "Alph_N_" + N.ToString() + "_Pop_" + NPop.ToString() + "_intensity" + intensity.ToString() + "alpha_" + alpha.ToString() + "_c_" + c.ToString() + method + "_Acc_" + Ac.Text + "JR" + Jr.ToString() + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_struct_" + structure + "OF_" + OF.ToString() + ".txt";
                                //if (!File.Exists(ss2))
                                //{
                                //    for (int tt = 0; tt < Ntest; tt++)
                                //    {
                                //        IntenI = new double[NPop];
                                //        AlphaI = new double[NPop];
                                //        JrI = new double[NPop];
                                //        IntenAr = new ArrayList();
                                //        AlphaAR = new ArrayList();
                                //        JrAR = new ArrayList();
                                //        NF[] Nodes = new NF[NPop];
                                //        NF[] ONodes = new NF[NPop];
                                //        NF[] Forces = new NF[NPop];
                                //        NF[] a = new NF[NPop];
                                //        NF[] v = new NF[NPop];

                                //        intensity = 0.5;
                                //        alpha = 0.5;
                                //        Jr = 0.3;
                                //        double[] Magnet = new double[NPop];
                                //        double[] Mass = new double[NPop];
                                //        NF Best_indi = new NF();
                                //        double[] distances = new double[NPop];
                                //        Best_indi.F = 1E290;
                                //        for (int ii = 0; ii < NPop; ii++)
                                //        {
                                //            Forces[ii].x = new double[N];
                                //            a[ii].x = new double[N];
                                //            v[ii].x = new double[N];
                                //        }
                                //        for (int ii = 0; ii < NPop; ii++)
                                //        {
                                //            Nodes[ii].x = NFC.node_constructing2(N, Rand,ub,lb);
                                //            Nodes[ii].F = class_new.Choosing_Function(OF, Nodes[ii].x, dimension, xopt, Perm, r25, r50, r100, S_Array, W_Array, overlap);

                                //            //PervFitness[ii] = Nodes[ii].F;
                                //            //if (OF < 16 || OF == 22)
                                //            //    ONodes[ii].F = NFC.compute_fitness(ONodes[ii].x, OF, N);
                                //            //else
                                //            //    ONodes[ii].F = NFC.compute_fitness(ONodes[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                //            //if (Nodes[ii].F < ONodes[ii].F)
                                //            //{
                                //            //    Nodes[ii].x = (double[])ONodes[ii].x.Clone();
                                //            //    Nodes[ii].F = ONodes[ii].F;
                                //            //}
                                //        }
                                //        for (int it = 0; it < Iteration; it++)
                                //        {
                                //            IntenAr.Clear(); JrAR.Clear(); AlphaAR.Clear();
                                //            Fmax = -1E290; Fmin = 1E290;
                                //            int min_index = 0; int max_index = 0;
                                //            for (int ii = 0; ii < NPop; ii++)
                                //            {
                                //                JrI[ii] = NFC.Normal_Distribution(Rand, Jr, 0.01);
                                //                PervFitness[ii] = Nodes[ii].F;
                                //                Nodes[ii].F = class_new.Choosing_Function(OF, Nodes[ii].x, dimension, xopt, Perm, r25, r50, r100, S_Array, W_Array, overlap);

                                //                if (it > 0)
                                //                    if (Nodes[ii].F < PervFitness[ii])
                                //                    {
                                //                        IntenAr.Add(IntenI[ii]);
                                //                        AlphaAR.Add(AlphaI[ii]);
                                //                    }
                                //                if (JrI[ii] > Rand.NextDouble())
                                //                {

                                //                    ONodes[ii].x = NFC.generating_jumping(Nodes, NPop, N, ii);
                                //                    ONodes[ii].F = class_new.Choosing_Function(OF, ONodes[ii].x, dimension, xopt, Perm, r25, r50, r100, S_Array, W_Array, overlap);
                                //                    if (Nodes[ii].F > ONodes[ii].F)
                                //                    {
                                //                        Nodes[ii].x = (double[])ONodes[ii].x.Clone();
                                //                        Nodes[ii].F = ONodes[ii].F;
                                //                        JrAR.Add(JrI[ii]);
                                //                    }
                                //                    Nodes[ii].x = NFC.Node_treatment(Nodes[ii].x, N, Rand);
                                //                }

                                //                if (Fmax < Nodes[ii].F)
                                //                {
                                //                    Fmax = Nodes[ii].F;
                                //                    max_index = ii;
                                //                }
                                //                if (Fmin > Nodes[ii].F)
                                //                {
                                //                    Fmin = Nodes[ii].F;
                                //                    min_index = ii;
                                //                }
                                //            }

                                //            if (Best_indi.F > Fmin)
                                //            {
                                //                Best_indi.F = Fmin;
                                //                Best_indi.x = (double[])Nodes[min_index].x.Clone();
                                //            }
                                //            double range = Fmax - Fmin;
                                //            if (range == 0)
                                //            {
                                //                Fmax = -1000000000; Fmin = 1000000000;
                                //                min_index = 0; max_index = 0;
                                //                for (int ii = 0; ii < NPop; ii++)
                                //                {
                                //                    Nodes[ii].x = NFC.node_constructing2(N, Rand,ub,lb);
                                //                    ONodes[ii].x = NFC.Opposition_Learning(Nodes[ii].x);
                                //                }
                                //            }
                                //            for (int ii = 0; ii < NPop; ii++)
                                //            {
                                //                IntenI[ii] = NFC.Normal_Distribution(Rand, intensity, 0.1);
                                //                AlphaI[ii] = NFC.Normal_Distribution(Rand, alpha, 0.1);
                                //                Magnet[ii] = (Fmax-Nodes[ii].F) / (range);
                                //                Mass[ii] = Magnet[ii] * IntenI[ii] + AlphaI[ii];
                                //            }
                                //            for (int ii = 0; ii < NPop; ii++)
                                //            {
                                //                Forces[ii].x = new double[N];
                                //            }
                                //            for (int ii = 0; ii < NPop; ii++)
                                //            {
                                //                int[] Neigbour = NFC.Structures(ii, NPop, structure);
                                //                for (int jj = 0; jj < Neigbour.Length; jj++)
                                //                {
                                //                    if (ii != Neigbour[jj])
                                //                    {
                                //                        if (D.Text == "1")
                                //                            distances[ii] = distance1(Nodes[ii].x, Nodes[Neigbour[jj]].x, N);
                                //                        //else
                                //                        //    distances[ii] = distance4(LHC[ii].LH, LHC[Neigbour[jj]].LH, N, P);
                                //                        if (distances[ii] > short_range)
                                //                        {
                                //                            for (int n = 0; n < N; n++)
                                //                                Forces[ii].x[n] = Forces[ii].x[n] + (Nodes[Neigbour[jj]].x[n] - Nodes[ii].x[n]) * Magnet[Neigbour[jj]] / distances[ii];

                                //                        }

                                //                    }
                                //                }
                                //                //for (int n = 0; n < N; n++)
                                //                //    for (int p = 0; p < P; p++)
                                //                //    {
                                //                //        Forces[ii].x[n, p] = Forces[ii].x[n, p] / (NPop - 1);
                                //                //        //Forces[ii].x[n, p] = Forces[ii].x[n, p] * rr.NextDouble();
                                //                //    }

                                //            }
                                //            if (acceleration == 0)
                                //            {
                                //                for (int ii = 0; ii < NPop; ii++)
                                //                {
                                //                    for (int n = 0; n < N; n++)
                                //                    {
                                //                        v[ii].x[n] = (Forces[ii].x[n] / Mass[ii]);
                                //                        Nodes[ii].x[n] = Nodes[ii].x[n] + v[ii].x[n];
                                //                    }

                                //                    Nodes[ii].x = NFC.Node_treatment(Nodes[ii].x, N, Rand);
                                //                }
                                //            }
                                //            else if (acceleration == 1)
                                //            {
                                //                for (int ii = 0; ii < NPop; ii++)
                                //                {
                                //                    for (int n = 0; n < N; n++)
                                //                    {
                                //                        a[ii].x[n] = (Forces[ii].x[n] / Mass[ii]);
                                //                        v[ii].x[n] = v[ii].x[n] * Rand.NextDouble() + a[ii].x[n];
                                //                        Nodes[ii].x[n] = Nodes[ii].x[n] + v[ii].x[n];
                                //                    }

                                //                    Nodes[ii].x = NFC.Node_treatment(Nodes[ii].x, N, Rand);
                                //                }
                                //            }
                                //            //////update///////////////////////////////////////
                                //            if (JrAR.Count > 0)
                                //            {
                                //                Jr = (1 - c) * Jr + c * NFC.mean_simple(JrAR);
                                //            }
                                //            if (it > 0)
                                //            {
                                //                if (AlphaAR.Count > 0)
                                //                {
                                //                    alpha = (1 - c) * alpha + c * NFC.mean_Lehmer(AlphaAR);
                                //                    intensity = (1 - c) * intensity + c * NFC.mean_Lehmer(IntenAr);
                                //                }
                                //            }
                                //            ///////////////////////////////////////////////////
                                //            Js[tt, it] = Jr;
                                //            Alphas[tt, it] = alpha;
                                //            Rhos[tt, it] = intensity;
                                //            fitnesses[it] = Best_indi.F;
                                //        }
                                //        costs[tt] = Best_indi.F;
                                //        //Nodes_test[tt].x = (double[])Best_indi.x.Clone();
                                //    }
                                //    //sw = new StreamWriter(ss2);
                                //    //for (int tt = 0; tt < Ntest; tt++)
                                //    //{
                                //    //    sw.Write(costs[tt].ToString() + " ");
                                //    //}
                                //    //sw.Close();
                                //    sw = new StreamWriter(ss1);
                                //    for (int tt = 0; tt < Ntest; tt++)
                                //    {
                                //        for (int it = 0; it < Iteration; it++)
                                //        {
                                //            sw.Write(Js[tt, it].ToString() + " ");
                                //        }
                                //        sw.WriteLine();
                                //    }
                                //    sw.Close();
                                //    sw = new StreamWriter(ss3);
                                //    for (int tt = 0; tt < Ntest; tt++)
                                //    {
                                //        for (int it = 0; it < Iteration; it++)
                                //        {
                                //            sw.Write(Rhos[tt, it].ToString() + " ");
                                //        }
                                //        sw.WriteLine();
                                //    }
                                //    sw.Close();
                                //    sw = new StreamWriter(ss4);
                                //    for (int tt = 0; tt < Ntest; tt++)
                                //    {
                                //        for (int it = 0; it < Iteration; it++)
                                //        {
                                //            sw.Write(Alphas[tt, it].ToString() + " ");
                                //        }
                                //        sw.WriteLine();
                                //    }
                                //    sw.Close();
                                //}
                            }
                            #endregion
                            #region FRFMOA
                            else if (method == "FRFMOA")
                            {
                                //////////////////////////////////Initial parameter values/////////////////
                                int acceleration = Convert.ToInt32(Ac.Text);
                                int Distance = Convert.ToInt32(D.Text);
                                double al1 = Convert.ToDouble(Al.Text);
                                double al2 = Convert.ToDouble(Al2.Text);
                                double al_step = Convert.ToDouble(Al_step.Text);
                                double in1 = Convert.ToDouble(Inten.Text);
                                double in2 = Convert.ToDouble(Inten2.Text);
                                double in_step = Convert.ToDouble(Inten_step.Text);
                                double short_range = 0.01;
                                int instance = 21;
                                int Nconf = 0;
                                int Npar = 3;
                                int Bwhole = Convert.ToInt32(Bud.Text);
                                int Bused = 0;
                                int Niter = Convert.ToInt32(Math.Floor(2 + Math.Log(Npar)));
                                int[] Bj = new int[Niter + 1];
                                Bj[0] = Bwhole / Niter;
                                Nconf = Bj[0] / 5;
                                Random RR = new Random();
                                float[,] Configuration = new float[Nconf, Npar];
                                double[] Cost = new double[Nconf];
                                R[,] Rha = new R[instance, Nconf];
                                for (int jj = 0; jj < Nconf; jj++)
                                {
                                    for (int kk = 0; kk < Npar - 1; kk++)
                                    {
                                        Configuration[jj, kk] = Convert.ToSingle(RR.NextDouble() * 100);
                                    }
                                    Configuration[jj, Npar - 1] = Convert.ToSingle(RR.NextDouble());
                                }
                                ///////////////////////////////////////////////////////
                                bool[] dis = new bool[Nconf]; // The ones have been discarded...
                                int dis_conf = 0;
                                int best_cand = -1;
                                int ins_count = 0;
                                double Mean = 0;
                                string ss = Directory.GetCurrentDirectory() + "\\FRFMOA.txt";
                                if (!File.Exists(ss))
                                {
                                    for (int IT = 0; IT < Niter; IT++)
                                    {
                                        Rha = new R[instance, Nconf];
                                        dis = new bool[Nconf];
                                        double[] orders = new double[Nconf];
                                        for (int II = 1; II < instance && dis_conf <= Nconf; II++)
                                        {
                                            OF = II;
                                            if (OF != 16)
                                            {
                                                StreamReader SRR = new StreamReader("orthog" + M.ToString() + ".txt");
                                                string all = SRR.ReadToEnd();
                                                string[] all_with_space = all.Split('\n', ',');
                                                Orthogonal_Matrix = new double[M, M];
                                                for (int ii = 0; ii < M; ii++)
                                                    for (int jj = 0; jj < M; jj++)
                                                    {
                                                        Orthogonal_Matrix[ii, jj] = Convert.ToDouble(all_with_space[ii + jj]);
                                                    }
                                            }
                                            else
                                            {
                                                StreamReader SRR = new StreamReader("orthog" + N.ToString() + ".txt");
                                                string all = SRR.ReadToEnd();
                                                string[] all_with_space = all.Split('\n', ',');
                                                Orthogonal_Matrix = new double[N, N];
                                                for (int ii = 0; ii < N; ii++)
                                                    for (int jj = 0; jj < N; jj++)
                                                    {
                                                        Orthogonal_Matrix[ii, jj] = Convert.ToDouble(all_with_space[ii + jj]);
                                                    }
                                            }
                                            if (OF != 14)
                                                for (int jjj = 0; jjj < Nconf; jjj++)
                                                {
                                                    if (!dis[jjj])
                                                    {
                                                        Mean = 0;
                                                        double intensity = Configuration[jjj, 0]; double alpha = Configuration[jjj, 1]; double Jr = Configuration[jjj, 2];
                                                        for (int tt = 0; tt < Ntest; tt++)
                                                        {
                                                            NF[] Nodes = new NF[NPop];
                                                            NF[] ONodes = new NF[NPop];
                                                            NF[] Forces = new NF[NPop];
                                                            NF[] a = new NF[NPop];
                                                            NF[] v = new NF[NPop];
                                                            double[] Magnet = new double[NPop];
                                                            double[] Mass = new double[NPop];
                                                            NF Best_indi = new NF();
                                                            double[] distances = new double[NPop];
                                                            Best_indi.F = -1E290;
                                                            for (int ii = 0; ii < NPop; ii++)
                                                            {
                                                                Forces[ii].x = new double[N];
                                                                a[ii].x = new double[N];
                                                                v[ii].x = new double[N];
                                                            }
                                                            for (int ii = 0; ii < NPop; ii++)
                                                            {
                                                                Nodes[ii].x = NFC.node_constructing(N, Rand);
                                                                if (OF < 16 || OF == 22)
                                                                    Nodes[ii].F = NFC.compute_fitness(Nodes[ii].x, OF, N);
                                                                else
                                                                    Nodes[ii].F = NFC.compute_fitness(Nodes[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                            }
                                                            for (int it = 0; it < Iteration; it++)
                                                            {
                                                                Fmax = -1E290; Fmin = 1000000000;
                                                                int min_index = 0; int max_index = 0;
                                                                for (int ii = 0; ii < NPop; ii++)
                                                                {
                                                                    if (OF < 16 || OF == 22)
                                                                        Nodes[ii].F = NFC.compute_fitness(Nodes[ii].x, OF, N);
                                                                    else
                                                                        Nodes[ii].F = NFC.compute_fitness(Nodes[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);

                                                                    if (Jr > Rand.NextDouble())
                                                                    {

                                                                        ONodes[ii].x = NFC.generating_jumping(Nodes, NPop, N, ii);
                                                                        if (OF < 16 || OF == 22)
                                                                            ONodes[ii].F = NFC.compute_fitness(ONodes[ii].x, OF, N);
                                                                        else
                                                                            ONodes[ii].F = NFC.compute_fitness(ONodes[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                                        if (Nodes[ii].F < ONodes[ii].F)
                                                                        {
                                                                            Nodes[ii].x = (double[])ONodes[ii].x.Clone();
                                                                            Nodes[ii].F = ONodes[ii].F;
                                                                        }
                                                                        Nodes[ii].x = NFC.Node_treatment(Nodes[ii].x, N, Rand);
                                                                    }

                                                                    if (Fmax < Nodes[ii].F)
                                                                    {
                                                                        Fmax = Nodes[ii].F;
                                                                        max_index = ii;
                                                                    }
                                                                    if (Fmin > Nodes[ii].F)
                                                                    {
                                                                        Fmin = Nodes[ii].F;
                                                                        min_index = ii;
                                                                    }
                                                                }

                                                                if (Best_indi.F < Fmax)
                                                                {
                                                                    Best_indi.F = Fmax;
                                                                    Best_indi.x = (double[])Nodes[max_index].x.Clone();
                                                                }
                                                                double range = Fmax - Fmin;
                                                                if (range == 0)
                                                                {
                                                                    Fmax = -1000000000; Fmin = 1000000000;
                                                                    min_index = 0; max_index = 0;
                                                                    for (int ii = 0; ii < NPop; ii++)
                                                                    {
                                                                        Nodes[ii].x = NFC.node_constructing(N, Rand);
                                                                        ONodes[ii].x = NFC.Opposition_Learning(Nodes[ii].x);
                                                                    }
                                                                }
                                                                for (int ii = 0; ii < NPop; ii++)
                                                                {
                                                                    Magnet[ii] = (Nodes[ii].F - Fmin) / (range);
                                                                    Mass[ii] = Magnet[ii] * intensity + alpha;
                                                                }
                                                                for (int ii = 0; ii < NPop; ii++)
                                                                {
                                                                    Forces[ii].x = new double[N];
                                                                }
                                                                for (int ii = 0; ii < NPop; ii++)
                                                                {
                                                                    int[] Neigbour = NFC.Structures(ii, NPop, structure);
                                                                    for (int jj = 0; jj < Neigbour.Length; jj++)
                                                                    {
                                                                        if (ii != Neigbour[jj])
                                                                        {
                                                                            if (D.Text == "1")
                                                                                distances[ii] = distance1(Nodes[ii].x, Nodes[Neigbour[jj]].x, N);
                                                                            //else
                                                                            //    distances[ii] = distance4(LHC[ii].LH, LHC[Neigbour[jj]].LH, N, P);
                                                                            if (distances[ii] > short_range)
                                                                            {
                                                                                for (int n = 0; n < N; n++)
                                                                                    Forces[ii].x[n] = Forces[ii].x[n] + (Nodes[Neigbour[jj]].x[n] - Nodes[ii].x[n]) * Magnet[Neigbour[jj]] / distances[ii];

                                                                            }

                                                                        }
                                                                    }

                                                                }
                                                                if (acceleration == 0)
                                                                {
                                                                    for (int ii = 0; ii < NPop; ii++)
                                                                    {
                                                                        for (int n = 0; n < N; n++)
                                                                        {
                                                                            v[ii].x[n] = (Forces[ii].x[n] / Mass[ii]);
                                                                            Nodes[ii].x[n] = Nodes[ii].x[n] + v[ii].x[n];
                                                                        }

                                                                        Nodes[ii].x = NFC.Node_treatment(Nodes[ii].x, N, Rand);
                                                                    }
                                                                }
                                                                else if (acceleration == 1)
                                                                {
                                                                    for (int ii = 0; ii < NPop; ii++)
                                                                    {
                                                                        for (int n = 0; n < N; n++)
                                                                        {
                                                                            a[ii].x[n] = (Forces[ii].x[n] / Mass[ii]);
                                                                            v[ii].x[n] = v[ii].x[n] * Rand.NextDouble() + a[ii].x[n];
                                                                            Nodes[ii].x[n] = Nodes[ii].x[n] + v[ii].x[n];
                                                                        }

                                                                        Nodes[ii].x = NFC.Node_treatment(Nodes[ii].x, N, Rand);
                                                                    }
                                                                }
                                                                //////update///////////////////////////////////////

                                                                ///////////////////////////////////////////////////

                                                            }
                                                            costs[tt] = Best_indi.F;
                                                            Mean += Best_indi.F;
                                                        }
                                                        Rha[II, jjj].Cost = Mean / Ntest;
                                                        Rha[II, jjj].conf = jjj;
                                                    }
                                                }
                                            for (int hh = II; hh < II + 1; hh++)
                                            {
                                                for (int hh2 = 0; hh2 < Nconf - 1; hh2++)
                                                {
                                                    if (!dis[hh2])
                                                        for (int hh3 = hh2 + 1; hh3 < Nconf; hh3++)
                                                        {
                                                            if (!dis[hh3])
                                                                if (Rha[hh, hh2].Cost < Rha[hh, hh3].Cost)
                                                                {
                                                                    R temp = Rha[hh, hh2];
                                                                    Rha[hh, hh2] = Rha[hh, hh3];
                                                                    Rha[hh, hh3] = temp;
                                                                }
                                                        }
                                                }
                                            }
                                            for (int hh = II; hh < II + 1; hh++)
                                            {
                                                for (int ii = Nconf - 1; ii > 0; ii--)
                                                {
                                                    if (!dis[ii])
                                                    {
                                                        Rha[hh, ii].value++;
                                                        for (int jj = ii - 1; jj >= 0; jj--)
                                                        {
                                                            if (Rha[hh, ii].Cost == Rha[hh, jj].Cost)
                                                            {
                                                                Rha[hh, ii].value++;
                                                            }
                                                            else
                                                            {
                                                                ii = jj + 1;
                                                                break;
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                            if (Rha[II, 1].value <= 1 && !dis[0])
                                            {
                                                Rha[II, 0].value = 1;
                                            }
                                            int Nc = 0;
                                            for (int ii = 0; ii < Nconf; ii++)
                                            {
                                                if (!dis[ii])
                                                {

                                                    if (Rha[II, ii].value > 2)
                                                    {
                                                        int temp = Convert.ToInt32(Rha[II, ii].value);
                                                        int ind = 0;
                                                        for (int jj = ii - temp + 1; jj <= ii; jj++)
                                                        {
                                                            ind += jj + 1;
                                                        }
                                                        for (int jj = ii - temp + 1; jj <= ii; jj++)
                                                        {
                                                            Rha[II, jj].value = ind / temp;
                                                        }
                                                    }
                                                    else
                                                    {
                                                        Rha[II, ii].value = ii + 1 - Nc;
                                                    }
                                                }
                                                else
                                                    Nc++;

                                            }
                                            for (int jj = II; jj < II + 1; jj++)
                                                for (int ii = 0; ii < orders.Length; ii++)
                                                {
                                                    if (!dis[ii])
                                                    {
                                                        orders[Rha[jj, ii].conf] += Rha[II, ii].value;
                                                    }
                                                }
                                            ////////////////////T computation
                                            double T = 0;
                                            double sorat = 0;
                                            for (int hh2 = 0; hh2 < Nconf; hh2++)
                                            {
                                                if (!dis[hh2])
                                                {
                                                    sorat += Math.Pow(orders[hh2] - (((II + 1) * Nconf + 1) / 2), 2);
                                                }
                                            }
                                            sorat = sorat * (Nconf - 1);
                                            double makhraj = 0;
                                            for (int hh = 0; hh < II + 1; hh++)
                                            {
                                                for (int hh2 = 0; hh2 < Nconf; hh2++)
                                                {
                                                    if (!dis[hh2])
                                                    {
                                                        makhraj += Math.Pow(Rha[hh, hh2].conf, 2);
                                                    }
                                                }
                                            }

                                            makhraj = makhraj - ((II * Nconf * Math.Pow((Nconf + 1), 2)) / 4);
                                            T = sorat / makhraj;
                                            ////////////////
                                            double alpha2 = 0.05;
                                            double quantile1 = NFC.get_t_quantile(alpha2, II + 1, t_table); //getStudentT(alpha, rows);
                                            if (T <= quantile1)
                                            {

                                            }
                                            else
                                            {
                                                //calculating T2 (pair-wise test), when needed
                                                double T2, quantile2 = NFC.get_t_quantile(alpha2, II + 1, t_table, true);
                                                double best_score = -1;
                                                for (int j = 0; j < Nconf; j++)
                                                {

                                                    if (orders[j] > best_score || best_score == -1)
                                                    {
                                                        best_score = orders[j];
                                                        best_cand = j;
                                                    }
                                                }
                                                makhraj *= Math.Abs((2 * (II + 1)) * (1 - (T / ((II + 1) * (Nconf - 1)))));
                                                makhraj /= ((II + 1) * (Nconf - 1));
                                                for (int j = 0; j < Nconf; j++)
                                                    if (!dis[j])
                                                    {
                                                        if (j != best_cand)
                                                        {
                                                            double diff = orders[best_cand] - orders[j];
                                                            if (diff < 0) diff *= (-1);
                                                            T2 = diff / Math.Sqrt(makhraj);
                                                            if (diff / Math.Sqrt(makhraj) > quantile2)
                                                            {
                                                                dis[j] = true;
                                                                dis_conf++;
                                                                orders[j] = 0;
                                                            }
                                                        }
                                                    }
                                            }
                                            ////////////////////
                                            ins_count = II;
                                        }
                                        //////////////////////////////////////////Select The best///////////////////////////////////
                                        int remainder = Nconf - dis_conf;
                                        int[] rz = new int[Nconf];
                                        int kk2 = 0;
                                        for (int ii = 0; ii < rz.Length; ii++)
                                        {
                                            if (!dis[ii])
                                            {
                                                rz[Rha[ins_count, ii].conf] = kk2 + 1;
                                                kk2++;
                                            }
                                        }
                                        double[] Pz = new double[remainder];
                                        double rands = RR.NextDouble();
                                        kk2 = 0;
                                        for (int ii = 0; ii < Nconf; ii++)
                                        {
                                            if (!dis[ii])
                                            {
                                                if (kk2 == 0)
                                                    Pz[kk2] = Convert.ToDouble(remainder - orders[ii] + 1) / (remainder * (remainder + 1) / 2);
                                                else
                                                    Pz[kk2] = Convert.ToDouble(remainder - orders[ii] + 1) / (remainder * (remainder + 1) / 2) + Pz[kk2 - 1];
                                                if (Pz[kk2] > rands)
                                                {
                                                    best_cand = Rha[ins_count, ii].conf;
                                                    break;
                                                }
                                                kk2++;
                                            }
                                        }
                                        /////////////////////////////////////////////////////////////////////////////////////////
                                        Bused += Bj[IT];
                                        double lboundCW = 0;
                                        double uboundCW = 10;
                                        double l_lg = 0;
                                        double u_lg = 1;
                                        int Mue = 5 + IT + 1;
                                        double deltad_menfi_yekCW = (uboundCW - lboundCW) / 2;
                                        double deltad_menfi_yek_lg = (u_lg - l_lg) / 2;
                                        int symbol = 0;
                                        Bj[IT + 1] = (Bwhole - Bused) / (Niter - IT + 1);
                                        int Nconf_perv = Nconf;
                                        Nconf = Convert.ToInt32(Bj[IT + 1] / (Mue + Math.Min(5, IT + 1))) - remainder;// (Nconf-dis_conf) is equal the remainder of perv_conf
                                        if (Nconf < 1)
                                        {
                                            break;
                                        }
                                        double deltad_CW = deltad_menfi_yekCW * Math.Pow((1.0 / Nconf), ((double)IT) / Npar);
                                        double delta_lg = deltad_menfi_yek_lg * Math.Pow((1.0 / Nconf), ((double)IT) / Npar);
                                        float[,] PervConf = (float[,])Configuration.Clone();
                                        Configuration = new float[Nconf + remainder, Npar];
                                        ///////C1 , C2 , W, r_lg
                                        for (int jj = 0; jj < Nconf; jj++)
                                        {
                                            for (kk2 = 0; kk2 < Npar - 1; kk2++)
                                            {
                                                symbol = RR.Next(0, 2);
                                                if (symbol == 0)
                                                    Configuration[jj, kk2] = Convert.ToSingle(PervConf[best_cand, kk2] + deltad_CW * RR.NextDouble());
                                                else
                                                    Configuration[jj, kk2] = Convert.ToSingle(PervConf[best_cand, kk2] - deltad_CW * RR.NextDouble());
                                            }
                                            symbol = RR.Next(0, 2);
                                            if (symbol == 0)
                                                Configuration[jj, Npar - 1] = Convert.ToSingle(PervConf[best_cand, Npar - 1] + delta_lg * RR.NextDouble());
                                            else
                                                Configuration[jj, Npar - 1] = Convert.ToSingle(PervConf[best_cand, Npar - 1] - delta_lg * RR.NextDouble());
                                        }
                                        int iii = Nconf;
                                        for (int jj = 0; jj < Nconf_perv; jj++)
                                        {
                                            if (!dis[jj])
                                            {
                                                for (kk2 = 0; kk2 < Npar; kk2++)
                                                {
                                                    Configuration[iii, kk2] = PervConf[jj, kk2];
                                                }
                                                iii++;
                                            }
                                        }
                                        Nconf = Nconf + remainder;// the set of perv_conf+ curr_conf
                                        dis_conf = 0;
                                    }
                                    sw = new StreamWriter(ss);
                                    for (int ii = 0; ii < Npar; ii++)
                                    {
                                        sw.WriteLine(Configuration[best_cand, ii]);
                                    }
                                    sw.Close();
                                }
                            }
                            #endregion
                            #region OMOA
                            if (method == "OMOA")
                            {
                                double Jr = Convert.ToDouble(Jumping_rate.Text);
                                double short_range = 0.01;
                                int acceleration = Convert.ToInt32(Ac.Text);
                                int Distance = Convert.ToInt32(D.Text);
                                double al1 = Convert.ToDouble(Al.Text);
                                double al2 = Convert.ToDouble(Al2.Text);
                                double al_step = Convert.ToDouble(Al_step.Text);
                                double in1 = Convert.ToDouble(Inten.Text);
                                double in2 = Convert.ToDouble(Inten2.Text);
                                double in_step = Convert.ToDouble(Inten_step.Text);
                                //double[] Intinsities = { 0.001, 0.01, 0.2, 0.6, 1 };
                                //double[] Alphas = { 0.001, 0.01, 0.2, 0.6, 1 };
                                //double[] Jrs = { 0.3, 0.45, 0.6 };
                                //double[] Intinsities = { 0.2, 0.001, 0.01, 0.2, 0.2, 0.01 };
                                //double[] Alphas = { 1, 0.6, 0.2, 0.001, 0.001, 0.2 };
                                //double[] Jrs = { 0.3, 0.45, 0.6, 0.6, 0.6, 0.6 };
                                //int ITT = 0;
                                double[,] fits = new double[Ntest, Iteration];
                                double[] Intinsities = { 0.6, 1.0, 1.0, 0.2, 0.6, 1.0, 0.2, 1.0, 1.0, 1.0, 1.0, 0.6, 1.0, 0, 1.0, 0.2, 0.2, 0.0010, 0.0010, 0.00100, 0.0010, 1.0 };
                                double[] Alphas = { 0.0200000000000000, 0.400000000000000, 0.400000000000000, 0.200000000000000, 1, 0.800000000000000, 0.200000000000000, 0.200000000000000, 0.6, 0.4, 0.4, 1, 0.2, 0, 0.6, 1, 0.6, 0.4, 0.0010, 0.01000, 0.2000, 1 };
                                double[] Jrs = { 0.6, 0.3, 0.45, 0.3, 0.3, 0.3, 0.45, 0.3, 0.3, 0.3, 0.3, 0.3, 0.45, 0, 0.3000, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.3 };
                                //int[] POP = { 5, 10, 25, 50, 100, 200, 250, 500 };
                                //int[] ITs = { 10000, 5000, 2000, 1000, 500, 250, 200, 100 };
                                //for (int al_index = 0; al_index < Alphas.Length; al_index++)
                                //    for (int int_index = 0; int_index < Intinsities.Length; int_index++)
                                //        for (int II = 0; II < 3; II++)
                                //        {
                                //double intensity = 93.01561;
                                //double alpha = 86.6102;
                                //Jr = 0.7958435;
                                double intensity = Intinsities[OF];
                                double alpha = Alphas[OF];
                                Jr = Jrs[OF];
                                //NPop = POP[II];
                                //Iteration = ITs[II];
                                //ss = "Result_N_" + N.ToString() + "_Pop_" + NPop.ToString() + "_intensity" + intensity.ToString() + "_alpha_" + alpha.ToString() + "_Dis_" + method + Distance.ToString() + "_Acc_" + Ac.Text + "JR" + Jr.ToString() + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_struct_" + structure + "OF_" + OF.ToString() + ".txt";
                                //ss2 = "Results\\fit_N_" + N.ToString() + "_bound_" + boundary_method.ToString() + "_Pop_" + NPop.ToString() + "_intensity" + intensity.ToString() + "alpha_" + alpha.ToString() + "_Dis_" + method + Distance.ToString() + "_Acc_" + Ac.Text + "JR" + Jr.ToString() + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_struct_" + structure + "OF_" + OF.ToString() + ".txt";
                                //string ss3 = "fits_N_" + N.ToString() + "_Pop_" + NPop.ToString() + "_intensity" + intensity.ToString() + "alpha_" + alpha.ToString() + "_Dis_" + method + Distance.ToString() + "_Acc_" + Ac.Text + "JR" + Jr.ToString() + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_struct_" + structure + "OF_" + OF.ToString() + ".txt";
                                if (!File.Exists(result))
                                {
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        NF[] Nodes = new NF[NPop];
                                        NF[] ONodes = new NF[NPop];
                                        NF[] Forces = new NF[NPop];
                                        NF[] a = new NF[NPop];
                                        NF[] v = new NF[NPop];
                                        double[] Magnet = new double[NPop];
                                        double[] Mass = new double[NPop];
                                        NF Best_indi = new NF();
                                        double[] distances = new double[NPop];
                                        Best_indi.F = -1E290;
                                        for (int ii = 0; ii < NPop; ii++)
                                        {
                                            Forces[ii].x = new double[N];
                                            a[ii].x = new double[N];
                                            v[ii].x = new double[N];
                                        }
                                        for (int ii = 0; ii < NPop; ii++)
                                        {
                                            Nodes[ii].x = NFC.node_constructing(N, Rand);
                                            ONodes[ii].x = NFC.Opposition_Learning(Nodes[ii].x);
                                        }
                                        for (int it = 0; it < Iteration; it++)
                                        {
                                            Fmax = -1E290; Fmin = 1000000000;
                                            int min_index = 0; int max_index = 0;
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                if (OF < 16 || OF == 22)
                                                    Nodes[ii].F = NFC.compute_fitness(Nodes[ii].x, OF, N);
                                                else
                                                    Nodes[ii].F = NFC.compute_fitness(Nodes[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);

                                                if (Jr > Rand.NextDouble())
                                                {

                                                    ONodes[ii].x = NFC.generating_jumping(Nodes, NPop, N, ii);
                                                    if (OF < 16 || OF == 22)
                                                        ONodes[ii].F = NFC.compute_fitness(ONodes[ii].x, OF, N);
                                                    else
                                                        ONodes[ii].F = NFC.compute_fitness(ONodes[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                    if (Nodes[ii].F < ONodes[ii].F)
                                                    {
                                                        Nodes[ii].x = (double[])ONodes[ii].x.Clone();
                                                        Nodes[ii].F = ONodes[ii].F;
                                                    }
                                                    Nodes[ii].x = NFC.Node_treatment(Nodes[ii].x, N, Rand);
                                                }

                                                if (Fmax < Nodes[ii].F)
                                                {
                                                    Fmax = Nodes[ii].F;
                                                    max_index = ii;
                                                }
                                                if (Fmin > Nodes[ii].F)
                                                {
                                                    Fmin = Nodes[ii].F;
                                                    min_index = ii;
                                                }
                                            }

                                            if (Best_indi.F < Fmax)
                                            {
                                                Best_indi.F = Fmax;
                                                Best_indi.x = (double[])Nodes[max_index].x.Clone();
                                            }
                                            double range = Fmax - Fmin;
                                            if (range == 0)
                                            {
                                                Fmax = -1000000000; Fmin = 1000000000;
                                                min_index = 0; max_index = 0;
                                                for (int ii = 0; ii < NPop; ii++)
                                                {
                                                    Nodes[ii].x = NFC.node_constructing(N, Rand);
                                                    ONodes[ii].x = NFC.Opposition_Learning(Nodes[ii].x);
                                                }
                                            }
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                Magnet[ii] = (Nodes[ii].F - Fmin) / (range);
                                                Mass[ii] = Magnet[ii] * intensity + alpha;
                                            }
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                Forces[ii].x = new double[N];
                                            }
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                int[] Neigbour = NFC.Structures(ii, NPop, structure);
                                                for (int jj = 0; jj < Neigbour.Length; jj++)
                                                {
                                                    if (ii != Neigbour[jj])
                                                    {
                                                        if (D.Text == "1")
                                                            distances[ii] = distance1(Nodes[ii].x, Nodes[Neigbour[jj]].x, N);
                                                        //else
                                                        //    distances[ii] = distance4(LHC[ii].LH, LHC[Neigbour[jj]].LH, N, P);
                                                        if (distances[ii] > short_range)
                                                        {
                                                            for (int n = 0; n < N; n++)
                                                                Forces[ii].x[n] = Forces[ii].x[n] + (Nodes[Neigbour[jj]].x[n] - Nodes[ii].x[n]) * Magnet[Neigbour[jj]] / distances[ii];

                                                        }

                                                    }
                                                }
                                                //for (int n = 0; n < N; n++)
                                                //    for (int p = 0; p < P; p++)
                                                //    {
                                                //        Forces[ii].x[n, p] = Forces[ii].x[n, p] / (NPop - 1);
                                                //        //Forces[ii].x[n, p] = Forces[ii].x[n, p] * rr.NextDouble();
                                                //    }

                                            }
                                            if (acceleration == 0)
                                            {
                                                for (int ii = 0; ii < NPop; ii++)
                                                {
                                                    for (int n = 0; n < N; n++)
                                                    {
                                                        v[ii].x[n] = (Forces[ii].x[n] / Mass[ii]);
                                                        Nodes[ii].x[n] = Nodes[ii].x[n] + v[ii].x[n];
                                                    }

                                                    //Nodes[ii].x = NFC.Node_treatment(Nodes[ii].x, N, Rand);
                                                    Solution_Treatment(boundary_method, ref Nodes[ii].x, Rand);
                                                }
                                            }
                                            else if (acceleration == 1)
                                            {
                                                for (int ii = 0; ii < NPop; ii++)
                                                {
                                                    for (int n = 0; n < N; n++)
                                                    {
                                                        a[ii].x[n] = (Forces[ii].x[n] / Mass[ii]);
                                                        v[ii].x[n] = v[ii].x[n] * Rand.NextDouble() + a[ii].x[n];
                                                        Nodes[ii].x[n] = Nodes[ii].x[n] + v[ii].x[n];
                                                    }

                                                    Nodes[ii].x = NFC.Node_treatment(Nodes[ii].x, N, Rand);
                                                }
                                            }

                                            fits[tt, it] = Best_indi.F;
                                        }
                                        costs[tt] = Best_indi.F;
                                        Nodes_test[tt].x = (double[])Best_indi.x.Clone();
                                    }
                                    sw = new StreamWriter(result);
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        sw.Write(costs[tt].ToString() + " ");
                                    }
                                    sw.Close();
                                    //sw = new StreamWriter(ss3);
                                    //for (int tt = 0; tt < Ntest; tt++)
                                    //{
                                    //    for (int ii = 0; ii < Iteration; ii++)
                                    //    {
                                    //        sw.Write(fits[tt, ii].ToString() + " ");
                                    //    }
                                    //    sw.WriteLine();
                                    //}
                                    //sw.Close();
                                }

                                //}
                            }
                            #endregion
                            #region FMOA
                            if (method == "FMOA")
                            {
                                //////////////////////////////////Initial parameter values/////////////////
                                double Jr = Convert.ToDouble(Jumping_rate.Text);
                                int acceleration = Convert.ToInt32(Ac.Text);
                                int Distance = Convert.ToInt32(D.Text);
                                double al1 = Convert.ToDouble(Al.Text);
                                double al2 = Convert.ToDouble(Al2.Text);
                                double al_step = Convert.ToDouble(Al_step.Text);
                                double in1 = Convert.ToDouble(Inten.Text);
                                double in2 = Convert.ToDouble(Inten2.Text);
                                double in_step = Convert.ToDouble(Inten_step.Text);
                                double short_range = 0.01;
                                //double[] Intinsities = { 0.001,0.01,0.2,0.6,1 };
                                //double[] Alphas = {0.001, 0.01,0.2,0.6,1 };
                                //double[] Jrs = { 0.3,0.45,0.6 };
                                //int ITT = 0;
                                //double[] Intinsities = { 0.6, 1.0, 1.0, 0.2, 0.6, 1.0, 0.2, 1.0, 1.0, 1.0, 1.0, 0.6, 1.0, 0, 1.0, 0.2, 0.2, 0.0010, 0.0010, 0.00100, 0.0010, 1.0 };
                                //double[] Alphas = { 0.0200000000000000, 0.400000000000000, 0.400000000000000, 0.200000000000000, 1, 0.800000000000000, 0.200000000000000, 0.200000000000000, 0.6, 0.4, 0.4, 1, 0.2, 0, 0.6, 1, 0.6, 0.4, 0.0010, 0.01000, 0.2000, 1 };
                                //double[] Jrs = { 0.6, 0.3, 0.45, 0.3, 0.3, 0.3, 0.45, 0.3, 0.3, 0.3, 0.3, 0.3, 0.45, 0, 0.3000, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.3 };
                                //int[] POP = { 5, 10, 25, 50, 100, 200, 250, 500 };
                                //int[] ITs = { 10000, 5000, 2000, 1000, 500, 250, 200, 100 };
                                //for (int al_index = 0; al_index < Alphas.Length; al_index++)
                                //    for (int int_index = 0; int_index < Intinsities.Length; int_index++)
                                //        for (int II = 0; II < 3; II++)
                                //        {
                                double intensity = 0.5;
                                double alpha = 0.5;
                                Jr = 0.3;
                                int maxall = 50000;
                                double c = 0.5;
                                double[] IntenI = new double[NPop];
                                double[] AlphaI = new double[NPop];
                                double[] JrI = new double[NPop];
                                ArrayList IntenAr = new ArrayList();
                                ArrayList AlphaAR = new ArrayList();
                                ArrayList JrAR = new ArrayList();
                                double[] PervFitness = new double[NPop];
                                double[,] Js = new double[Ntest, Iteration];
                                double[,] Rhos = new double[Ntest, Iteration];
                                double[,] Alphas = new double[Ntest, Iteration];
                                double[,] fitnesses = new double[Ntest, Iteration];
                                //NPop = POP[II];
                                //Iteration = ITs[II];
                                //for (int np = 0; np < Pops.Length; np++)
                                //{
                                    intensity = 0.5;
                                    alpha = 0.5;
                                    Jr = 0.3;
                                    //NPop = Pops[np];
                                    Iteration = maxall / NPop;
                                    //ss2 = "Results\\fit_N_" + N.ToString() + "_bound_" + boundary_method.ToString() +"_Pop_" + NPop.ToString() + "_intensity" + intensity.ToString() + "alpha_" + alpha.ToString() + "_c_" + c.ToString() + method + "_Acc_" + Ac.Text + "JR" + Jr.ToString() + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_struct_" + structure + "OF_" + OF.ToString() + ".txt";
                                    //string ss1 = "JR_N_" + N.ToString() + "_Pop_" + NPop.ToString() + "_intensity" + intensity.ToString() + "alpha_" + alpha.ToString() + "_c_" + c.ToString() + method + "_Acc_" + Ac.Text + "JR" + Jr.ToString() + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_struct_" + structure + "OF_" + OF.ToString() + ".txt";
                                    //string ss3 = "Rho_N_" + N.ToString() + "_Pop_" + NPop.ToString() + "_intensity" + intensity.ToString() + "alpha_" + alpha.ToString() + "_c_" + c.ToString() + method + "_Acc_" + Ac.Text + "JR" + Jr.ToString() + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_struct_" + structure + "OF_" + OF.ToString() + ".txt";
                                    //string ss4 = "Alph_N_" + N.ToString() + "_Pop_" + NPop.ToString() + "_intensity" + intensity.ToString() + "alpha_" + alpha.ToString() + "_c_" + c.ToString() + method + "_Acc_" + Ac.Text + "JR" + Jr.ToString() + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_struct_" + structure + "OF_" + OF.ToString() + ".txt";
                                    if (!File.Exists(result))
                                    {
                                        for (int tt = 0; tt < Ntest; tt++)
                                        {
                                            IntenI = new double[NPop];
                                            AlphaI = new double[NPop];
                                            JrI = new double[NPop];
                                            IntenAr = new ArrayList();
                                            AlphaAR = new ArrayList();
                                            JrAR = new ArrayList();
                                            NF[] Nodes = new NF[NPop];
                                            NF[] ONodes = new NF[NPop];
                                            NF[] Forces = new NF[NPop];
                                            NF[] a = new NF[NPop];
                                            NF[] v = new NF[NPop];

                                            intensity = 0.5;
                                            alpha = 0.5;
                                            Jr = 0.3;
                                            double[] Magnet = new double[NPop];
                                            double[] Mass = new double[NPop];
                                            NF Best_indi = new NF();
                                            double[] distances = new double[NPop];
                                            Best_indi.F = -1E290;
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                Forces[ii].x = new double[N];
                                                a[ii].x = new double[N];
                                                v[ii].x = new double[N];
                                            }
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                Nodes[ii].x = NFC.node_constructing(N, Rand);
                                                //ONodes[ii].x = NFC.Opposition_Learning(Nodes[ii].x);
                                                if (OF < 16 || OF == 22)
                                                    Nodes[ii].F = NFC.compute_fitness(Nodes[ii].x, OF, N);
                                                else
                                                    Nodes[ii].F = NFC.compute_fitness(Nodes[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                //PervFitness[ii] = Nodes[ii].F;
                                                //if (OF < 16 || OF == 22)
                                                //    ONodes[ii].F = NFC.compute_fitness(ONodes[ii].x, OF, N);
                                                //else
                                                //    ONodes[ii].F = NFC.compute_fitness(ONodes[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                //if (Nodes[ii].F < ONodes[ii].F)
                                                //{
                                                //    Nodes[ii].x = (double[])ONodes[ii].x.Clone();
                                                //    Nodes[ii].F = ONodes[ii].F;
                                                //}
                                            }
                                            for (int it = 0; it < Iteration; it++)
                                            {
                                                IntenAr.Clear(); JrAR.Clear(); AlphaAR.Clear();
                                                Fmax = -1E290; Fmin = 1000000000;
                                                int min_index = 0; int max_index = 0;
                                                for (int ii = 0; ii < NPop; ii++)
                                                {
                                                    JrI[ii] = NFC.Normal_Distribution(Rand, Jr, 0.01);
                                                    PervFitness[ii] = Nodes[ii].F;
                                                    if (OF < 16 || OF == 22)
                                                        Nodes[ii].F = NFC.compute_fitness(Nodes[ii].x, OF, N);
                                                    else
                                                        Nodes[ii].F = NFC.compute_fitness(Nodes[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                    if (it > 0)
                                                        if (Nodes[ii].F > PervFitness[ii])
                                                        {
                                                            IntenAr.Add(IntenI[ii]);
                                                            AlphaAR.Add(AlphaI[ii]);
                                                        }
                                                    if (JrI[ii] > Rand.NextDouble())
                                                    {

                                                        ONodes[ii].x = NFC.generating_jumping(Nodes, NPop, N, ii);
                                                        if (OF < 16 || OF == 22)
                                                            ONodes[ii].F = NFC.compute_fitness(ONodes[ii].x, OF, N);
                                                        else
                                                            ONodes[ii].F = NFC.compute_fitness(ONodes[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                        if (Nodes[ii].F < ONodes[ii].F)
                                                        {
                                                            Nodes[ii].x = (double[])ONodes[ii].x.Clone();
                                                            Nodes[ii].F = ONodes[ii].F;
                                                            JrAR.Add(JrI[ii]);
                                                        }
                                                        //Nodes[ii].x = NFC.Node_treatment(Nodes[ii].x, N, Rand);
                                                        Solution_Treatment(boundary_method, ref Nodes[ii].x, Rand);
                                                    }

                                                    if (Fmax < Nodes[ii].F)
                                                    {
                                                        Fmax = Nodes[ii].F;
                                                        max_index = ii;
                                                    }
                                                    if (Fmin > Nodes[ii].F)
                                                    {
                                                        Fmin = Nodes[ii].F;
                                                        min_index = ii;
                                                    }
                                                }

                                                if (Best_indi.F < Fmax)
                                                {
                                                    Best_indi.F = Fmax;
                                                    Best_indi.x = (double[])Nodes[max_index].x.Clone();
                                                }
                                                double range = Fmax - Fmin;
                                                if (range == 0)
                                                {
                                                    Fmax = -1000000000; Fmin = 1000000000;
                                                    min_index = 0; max_index = 0;
                                                    for (int ii = 0; ii < NPop; ii++)
                                                    {
                                                        Nodes[ii].x = NFC.node_constructing(N, Rand);
                                                        ONodes[ii].x = NFC.Opposition_Learning(Nodes[ii].x);
                                                    }
                                                }
                                                for (int ii = 0; ii < NPop; ii++)
                                                {
                                                    IntenI[ii] = NFC.Normal_Distribution(Rand, intensity, 0.1);
                                                    AlphaI[ii] = NFC.Normal_Distribution(Rand, alpha, 0.1);
                                                    Magnet[ii] = (Nodes[ii].F - Fmin) / (range);
                                                    Mass[ii] = Magnet[ii] * IntenI[ii] + AlphaI[ii];
                                                }
                                                for (int ii = 0; ii < NPop; ii++)
                                                {
                                                    Forces[ii].x = new double[N];
                                                }
                                                for (int ii = 0; ii < NPop; ii++)
                                                {
                                                    int[] Neigbour = NFC.Structures(ii, NPop, structure);
                                                    for (int jj = 0; jj < Neigbour.Length; jj++)
                                                    {
                                                        if (ii != Neigbour[jj])
                                                        {
                                                            if (D.Text == "1")
                                                                distances[ii] = distance1(Nodes[ii].x, Nodes[Neigbour[jj]].x, N);
                                                            //else
                                                            //    distances[ii] = distance4(LHC[ii].LH, LHC[Neigbour[jj]].LH, N, P);
                                                            if (distances[ii] > short_range)
                                                            {
                                                                for (int n = 0; n < N; n++)
                                                                    Forces[ii].x[n] = Forces[ii].x[n] + (Nodes[Neigbour[jj]].x[n] - Nodes[ii].x[n]) * Magnet[Neigbour[jj]] / distances[ii];

                                                            }

                                                        }
                                                    }
                                                    //for (int n = 0; n < N; n++)
                                                    //    for (int p = 0; p < P; p++)
                                                    //    {
                                                    //        Forces[ii].x[n, p] = Forces[ii].x[n, p] / (NPop - 1);
                                                    //        //Forces[ii].x[n, p] = Forces[ii].x[n, p] * rr.NextDouble();
                                                    //    }

                                                }
                                                if (acceleration == 0)
                                                {
                                                    for (int ii = 0; ii < NPop; ii++)
                                                    {
                                                        for (int n = 0; n < N; n++)
                                                        {
                                                            v[ii].x[n] = (Forces[ii].x[n] / Mass[ii]);
                                                            Nodes[ii].x[n] = Nodes[ii].x[n] + v[ii].x[n];
                                                        }

                                                        Nodes[ii].x = NFC.Node_treatment(Nodes[ii].x, N, Rand);
                                                    }
                                                }
                                                else if (acceleration == 1)
                                                {
                                                    for (int ii = 0; ii < NPop; ii++)
                                                    {
                                                        for (int n = 0; n < N; n++)
                                                        {
                                                            a[ii].x[n] = (Forces[ii].x[n] / Mass[ii]);
                                                            v[ii].x[n] = v[ii].x[n] * Rand.NextDouble() + a[ii].x[n];
                                                            Nodes[ii].x[n] = Nodes[ii].x[n] + v[ii].x[n];
                                                        }

                                                        //Nodes[ii].x = NFC.Node_treatment(Nodes[ii].x, N, Rand);
                                                        Solution_Treatment(boundary_method, ref Nodes[ii].x, Rand);
                                                    }
                                                }
                                                //////update///////////////////////////////////////
                                                if (JrAR.Count > 0)
                                                {
                                                    Jr = (1 - c) * Jr + c * NFC.mean_simple(JrAR);
                                                }
                                                if (it > 0)
                                                {
                                                    if (AlphaAR.Count > 0)
                                                    {
                                                        alpha = (1 - c) * alpha + c * NFC.mean_Lehmer(AlphaAR);
                                                        intensity = (1 - c) * intensity + c * NFC.mean_Lehmer(IntenAr);
                                                    }
                                                }
                                                /////////////////////////////////////////////////////
                                                //Js[tt, it] = Jr;
                                                //Alphas[tt, it] = alpha;
                                                //Rhos[tt, it] = intensity;
                                                //fitnesses[tt, it] = Best_indi.F;
                                            }
                                            costs[tt] = Best_indi.F;
                                            //Nodes_test[tt].x = (double[])Best_indi.x.Clone();
                                        }
                                        sw = new StreamWriter(result);
                                        for (int tt = 0; tt < Ntest; tt++)
                                        {
                                            sw.Write(costs[tt].ToString() + " ");
                                        }
                                        sw.Close();
                                        //sw = new StreamWriter(ss1);
                                        //for (int tt = 0; tt < Ntest; tt++)
                                        //{
                                        //    for (int it = 0; it < Iteration; it++)
                                        //    {
                                        //        sw.Write(Js[tt, it].ToString() + " ");
                                        //    }
                                        //    sw.WriteLine();
                                        //}
                                        //sw.Close();
                                        //sw = new StreamWriter(ss3);
                                        //for (int tt = 0; tt < Ntest; tt++)
                                        //{
                                        //    for (int it = 0; it < Iteration; it++)
                                        //    {
                                        //        sw.Write(Rhos[tt, it].ToString() + " ");
                                        //    }
                                        //    sw.WriteLine();
                                        //}
                                        //sw.Close();
                                        //sw = new StreamWriter(ss2);
                                        //for (int tt = 0; tt < Ntest; tt++)
                                        //{
                                        //    for (int it = 0; it < Iteration; it++)
                                        //    {
                                        //        sw.Write(fitnesses[tt, it].ToString() + " ");
                                        //    }
                                        //    sw.WriteLine();
                                        //}
                                        //sw.Close();
                                    }
                                //}
                            }
                            #endregion
                            #region NFFMOA
                            else if (method == "NFMOA")
                            {
                                //function_evaluation2013 class_new = new function_evaluation2013(ub, lb, dimension);

                                ////////////////////////////////////Initial parameter values/////////////////
                                //double Jr = Convert.ToDouble(Jumping_rate.Text);
                                //int acceleration = Convert.ToInt32(Ac.Text);
                                //int Distance = Convert.ToInt32(D.Text);
                                //double al1 = Convert.ToDouble(Al.Text);
                                //double al2 = Convert.ToDouble(Al2.Text);
                                //double al_step = Convert.ToDouble(Al_step.Text);
                                //double in1 = Convert.ToDouble(Inten.Text);
                                //double in2 = Convert.ToDouble(Inten2.Text);
                                //double in_step = Convert.ToDouble(Inten_step.Text);
                                //double short_range = 0.01;
                                ////double[] Intinsities = { 0.001,0.01,0.2,0.6,1 };
                                ////double[] Alphas = {0.001, 0.01,0.2,0.6,1 };
                                ////double[] Jrs = { 0.3,0.45,0.6 };
                                ////int ITT = 0;
                                ////double[] Intinsities = { 0.6, 1.0, 1.0, 0.2, 0.6, 1.0, 0.2, 1.0, 1.0, 1.0, 1.0, 0.6, 1.0, 0, 1.0, 0.2, 0.2, 0.0010, 0.0010, 0.00100, 0.0010, 1.0 };
                                ////double[] Alphas = { 0.0200000000000000, 0.400000000000000, 0.400000000000000, 0.200000000000000, 1, 0.800000000000000, 0.200000000000000, 0.200000000000000, 0.6, 0.4, 0.4, 1, 0.2, 0, 0.6, 1, 0.6, 0.4, 0.0010, 0.01000, 0.2000, 1 };
                                ////double[] Jrs = { 0.6, 0.3, 0.45, 0.3, 0.3, 0.3, 0.45, 0.3, 0.3, 0.3, 0.3, 0.3, 0.45, 0, 0.3000, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.3 };
                                ////int[] POP = { 5, 10, 25, 50, 100, 200, 250, 500 };
                                ////int[] ITs = { 10000, 5000, 2000, 1000, 500, 250, 200, 100 };
                                ////for (int al_index = 0; al_index < Alphas.Length; al_index++)
                                ////    for (int int_index = 0; int_index < Intinsities.Length; int_index++)
                                ////        for (int II = 0; II < 3; II++)
                                ////        {

                                //double[] adaptive_par = new double[] { 0.00, 0.05, 0.01, 0.05, 0.1,0.2,0.3,0.4 };
                                //double[] holded_iteration = new double[] { 0.001,0.005,0.01, 0.05, 0.1, 0.2 };//the percentage of iterations that the algorithm does not change its adaptive parameter.
                                //double[] fitnesses = new double[Iteration];
                                //double intensity = 0.5;
                                //double alpha = 0.5;
                                //Jr = 0.3;
                                //double c = 0.5;
                                //double[] IntenI = new double[NPop];
                                //double[] AlphaI = new double[NPop];
                                //double[] JrI = new double[NPop];
                                ////ArrayList IntenAr = new ArrayList();
                                ////ArrayList AlphaAR = new ArrayList();
                                ////ArrayList JrAR = new ArrayList();
                                //double[] PervFitness = new double[NPop];
                                //double[] ranges = new double[] { 0.001, 0.005, 0.01, 0.05, 0.1 };
                                ////NPop = POP[II];
                                ////Iteration = ITs[II];
                                //for (int ji = 0; ji < holded_iteration.Length;ji++ )
                                //    for (int ij = 0; ij < adaptive_par.Length; ij++)
                                //    {
                                //        int number_iteration_kept_in =Convert.ToInt32(holded_iteration[ji] * Iteration);
                                //        //ss2 = "fit_N_" + N.ToString() + "_Pop_" + NPop.ToString() + "_intensity" + intensity.ToString() + "alpha_" + alpha.ToString() + "_c_" + c.ToString() + method + "_Acc_" + Ac.Text + "JR" + Jr.ToString() + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_struct_" + structure + "OF_" + OF.ToString() + ".txt";
                                //        string ss1 = "JR_N_" + N.ToString() + "_Pop_" + NPop.ToString() + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "OF_" + "adap_par" + adaptive_par[ij].ToString() + OF.ToString()+"holded_iter"+holded_iteration[ji].ToString() + ".txt";
                                //        //string ss3 = "Rho_N_" + N.ToString() + "_Pop_" + NPop.ToString() + "_intensity" + intensity.ToString() + "alpha_" + alpha.ToString() + "_c_" + c.ToString() + method + "_Acc_" + Ac.Text + "JR" + Jr.ToString() + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_struct_" + structure + "OF_" + OF.ToString() + ".txt";
                                //        //string ss4 = "Alph_N_" + N.ToString() + "_Pop_" + NPop.ToString() + "_intensity" + intensity.ToString() + "alpha_" + alpha.ToString() + "_c_" + c.ToString() + method + "_Acc_" + Ac.Text + "JR" + Jr.ToString() + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_struct_" + structure + "OF_" + OF.ToString() + ".txt";
                                //        if (!File.Exists(ss2))
                                //        {
                                //            for (int tt = 0; tt < Ntest; tt++)
                                //            {
                                //                IntenI = new double[NPop];
                                //                AlphaI = new double[NPop];
                                //                //JrI = new double[NPop];
                                //                //IntenAr = new ArrayList();
                                //                //AlphaAR = new ArrayList();
                                //                //JrAR = new ArrayList();
                                //                NF[] Nodes = new NF[NPop];
                                //                NF[] ONodes = new NF[NPop];
                                //                NF[] Forces = new NF[NPop];
                                //                NF[] a = new NF[NPop];
                                //                NF[] v = new NF[NPop];
                                //                intensity = 0.5;
                                //                alpha = 0.5;
                                //                Jr = 0.3;

                                //                double[] Magnet = new double[NPop];
                                //                double[] Mass = new double[NPop];
                                //                NF Best_indi = new NF();
                                //                double[] distances = new double[NPop];
                                //                Best_indi.F = 1E290;

                                //                for (int ii = 0; ii < NPop; ii++)
                                //                {
                                //                    Forces[ii].x = new double[N];
                                //                    a[ii].x = new double[N];
                                //                    v[ii].x = new double[N];
                                //                    JrI[ii] = Jr;
                                //                    IntenI[ii] = intensity;
                                //                    AlphaI[ii] = alpha;
                                //                }
                                //                for (int ii = 0; ii < NPop; ii++)
                                //                {
                                //                    Nodes[ii].x = NFC.node_constructing2(N, Rand, ub, lb);
                                //                    Nodes[ii].F = class_new.Choosing_Function(OF, Nodes[ii].x, dimension, xopt, Perm, r25, r50, r100, S_Array, W_Array, overlap);
                                //                }
                                //                for (int it = 0; it < Iteration; it++)
                                //                {
                                //                    Fmax = -1E290; Fmin = 1E290;
                                //                    int min_index = 0; int max_index = 0;
                                //                    for (int ii = 0; ii < NPop; ii++)
                                //                    {
                                //                        JrI[ii] = NFC.Normal_Distribution(Rand, Jr, 0.01);
                                //                        PervFitness[ii] = Nodes[ii].F;
                                //                        Nodes[ii].F = class_new.Choosing_Function(OF, Nodes[ii].x, dimension, xopt, Perm, r25, r50, r100, S_Array, W_Array, overlap);

                                //                        if (it > 0)
                                //                            if (Nodes[ii].F < PervFitness[ii])
                                //                            {
                                //                                if ((it % number_iteration_kept_in) == 0)
                                //                                {
                                //                                    IntenI[ii] = IntenI[ii] + adaptive_par[ij] * Rand.NextDouble();
                                //                                    AlphaI[ii] = AlphaI[ii] + adaptive_par[ij] * Rand.NextDouble();
                                //                                }
                                //                            }
                                //                            else if (Nodes[ii].F > PervFitness[ii])
                                //                            {
                                //                                if ((it % number_iteration_kept_in) == 0)
                                //                                {
                                //                                    IntenI[ii] = IntenI[ii] - adaptive_par[ij] * Rand.NextDouble();
                                //                                    AlphaI[ii] = AlphaI[ii] - adaptive_par[ij] * Rand.NextDouble();
                                //                                }
                                //                            }
                                //                        //if (JrI[ii] > Rand.NextDouble())
                                //                        //{

                                //                        //    ONodes[ii].x = NFC.generating_jumping(Nodes, NPop, N, ii);
                                //                        //    ONodes[ii].F = class_new.Choosing_Function(OF, ONodes[ii].x, dimension, xopt, Perm, r25, r50, r100, S_Array, W_Array, overlap);
                                //                        //    if (Nodes[ii].F > ONodes[ii].F)
                                //                        //    {
                                //                        //        Nodes[ii].x = (double[])ONodes[ii].x.Clone();
                                //                        //        Nodes[ii].F = ONodes[ii].F;
                                //                        //        JrI[ii] = JrI[ii] + 0.01 * Rand.NextDouble();
                                //                        //    }
                                //                        //    else if (Nodes[ii].F < ONodes[ii].F)
                                //                        //        JrI[ii] = JrI[ii] - 0.01 * Rand.NextDouble();
                                //                        //    Nodes[ii].x = NFC.Node_treatment(Nodes[ii].x, N, Rand);
                                //                        //}

                                //                        if (Fmax < Nodes[ii].F)
                                //                        {
                                //                            Fmax = Nodes[ii].F;
                                //                            max_index = ii;
                                //                        }
                                //                        if (Fmin > Nodes[ii].F)
                                //                        {
                                //                            Fmin = Nodes[ii].F;
                                //                            min_index = ii;
                                //                        }
                                //                    }

                                //                    if (Best_indi.F > Fmin)
                                //                    {
                                //                        Best_indi.F = Fmin;
                                //                        Best_indi.x = (double[])Nodes[min_index].x.Clone();
                                //                    }
                                //                    double range = Fmax - Fmin;
                                //                    if (range == 0)
                                //                    {
                                //                        Fmax = -1000000000; Fmin = 1000000000;
                                //                        min_index = 0; max_index = 0;
                                //                        for (int ii = 0; ii < NPop; ii++)
                                //                        {
                                //                            Nodes[ii].x = NFC.node_constructing2(N, Rand, ub, lb);
                                //                            ONodes[ii].x = NFC.Opposition_Learning(Nodes[ii].x);
                                //                        }
                                //                    }
                                //                    for (int ii = 0; ii < NPop; ii++)
                                //                    {
                                //                        Magnet[ii] = (Fmax - Nodes[ii].F) / (range);
                                //                        Mass[ii] = Magnet[ii] * IntenI[ii] + AlphaI[ii];
                                //                    }
                                //                    for (int ii = 0; ii < NPop; ii++)
                                //                    {
                                //                        Forces[ii].x = new double[N];
                                //                    }
                                //                    for (int ii = 0; ii < NPop; ii++)
                                //                    {
                                //                        int[] Neigbour = NFC.Structures(ii, NPop, structure);
                                //                        for (int jj = 0; jj < Neigbour.Length; jj++)
                                //                        {
                                //                            if (ii != Neigbour[jj])
                                //                            {
                                //                                if (D.Text == "1")
                                //                                    distances[ii] = distance1(Nodes[ii].x, Nodes[Neigbour[jj]].x, N);
                                //                                //else
                                //                                //    distances[ii] = distance4(LHC[ii].LH, LHC[Neigbour[jj]].LH, N, P);
                                //                                if (distances[ii] > short_range)
                                //                                {
                                //                                    for (int n = 0; n < N; n++)
                                //                                        Forces[ii].x[n] = Forces[ii].x[n] + (Nodes[Neigbour[jj]].x[n] - Nodes[ii].x[n]) * Magnet[Neigbour[jj]] / distances[ii];

                                //                                }

                                //                            }
                                //                        }

                                //                    }
                                //                    if (acceleration == 0)
                                //                    {
                                //                        for (int ii = 0; ii < NPop; ii++)
                                //                        {
                                //                            for (int n = 0; n < N; n++)
                                //                            {
                                //                                v[ii].x[n] = (Forces[ii].x[n] / Mass[ii]);
                                //                                Nodes[ii].x[n] = Nodes[ii].x[n] + v[ii].x[n];
                                //                            }

                                //                            Nodes[ii].x = NFC.Node_treatment(Nodes[ii].x, N, Rand);
                                //                        }
                                //                    }
                                //                    else if (acceleration == 1)
                                //                    {
                                //                        for (int ii = 0; ii < NPop; ii++)
                                //                        {
                                //                            for (int n = 0; n < N; n++)
                                //                            {
                                //                                a[ii].x[n] = (Forces[ii].x[n] / Mass[ii]);
                                //                                v[ii].x[n] = v[ii].x[n] * Rand.NextDouble() + a[ii].x[n];
                                //                                Nodes[ii].x[n] = Nodes[ii].x[n] + v[ii].x[n];
                                //                            }

                                //                            Nodes[ii].x = NFC.Node_treatment(Nodes[ii].x, N, Rand);
                                //                        }
                                //                    }
                                //                    ///////////////////////////////////////////////////

                                //                    fitnesses[it] = Best_indi.F;
                                //                }
                                //                costs[tt] = Best_indi.F;
                                //                //Nodes_test[tt].x = (double[])Best_indi.x.Clone();
                                //            }
                                //            sw = new StreamWriter(ss1);
                                //            for (int tt = 0; tt < Ntest; tt++)
                                //            {
                                //                    sw.Write(costs[tt].ToString() + " ");
                                //            }
                                //            sw.Close();
                                //        }
                                //    }
                            }
                            #endregion
                            #region FMOAT
                            if (method == "FMOAT")
                            {
                                //////////////////////////////////Initial parameter values/////////////////
                                //double Jr = Convert.ToDouble(Jumping_rate.Text);
                                //int acceleration = Convert.ToInt32(Ac.Text);
                                //int Distance = Convert.ToInt32(D.Text);
                                //double al1 = Convert.ToDouble(Al.Text);
                                //double al2 = Convert.ToDouble(Al2.Text);
                                //double al_step = Convert.ToDouble(Al_step.Text);
                                //double in1 = Convert.ToDouble(Inten.Text);
                                //double in2 = Convert.ToDouble(Inten2.Text);
                                //double in_step = Convert.ToDouble(Inten_step.Text);
                                //double short_range = 0.01;
                                ////double[] Intinsities = { 0.001,0.01,0.2,0.6,1 };
                                ////double[] Alphas = {0.001, 0.01,0.2,0.6,1 };
                                ////double[] Jrs = { 0.3,0.45,0.6 };
                                ////int ITT = 0;
                                ////double[] Intinsities = { 0.6, 1.0, 1.0, 0.2, 0.6, 1.0, 0.2, 1.0, 1.0, 1.0, 1.0, 0.6, 1.0, 0, 1.0, 0.2, 0.2, 0.0010, 0.0010, 0.00100, 0.0010, 1.0 };
                                ////double[] Alphas = { 0.0200000000000000, 0.400000000000000, 0.400000000000000, 0.200000000000000, 1, 0.800000000000000, 0.200000000000000, 0.200000000000000, 0.6, 0.4, 0.4, 1, 0.2, 0, 0.6, 1, 0.6, 0.4, 0.0010, 0.01000, 0.2000, 1 };
                                ////double[] Jrs = { 0.6, 0.3, 0.45, 0.3, 0.3, 0.3, 0.45, 0.3, 0.3, 0.3, 0.3, 0.3, 0.45, 0, 0.3000, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.3 };
                                ////int[] POP = { 5, 10, 25, 50, 100, 200, 250, 500 };
                                ////int[] ITs = { 10000, 5000, 2000, 1000, 500, 250, 200, 100 };
                                ////for (int al_index = 0; al_index < Alphas.Length; al_index++)
                                ////    for (int int_index = 0; int_index < Intinsities.Length; int_index++)
                                ////        for (int II = 0; II < 3; II++)
                                ////        {
                                //double intensity = 0.5;
                                //double alpha = 0.5;
                                //Jr = 0.0;
                                //double c = 0.5;
                                //double[] IntenI = new double[NPop];
                                //double[] AlphaI = new double[NPop];
                                //double[] JrI = new double[NPop];
                                //ArrayList IntenAr = new ArrayList();
                                //ArrayList AlphaAR = new ArrayList();
                                //ArrayList JrAR = new ArrayList();
                                //double[] PervFitness = new double[NPop];
                                //int iter = Iteration / 100;
                                //double[,] Js = new double[Ntest, Iteration];
                                //double[,] Rhos = new double[Ntest, Iteration];
                                //double[,] Alphas = new double[Ntest, Iteration];
                                //double[,] Fs = new double[Ntest, Iteration];
                                ////NPop = POP[II];
                                ////Iteration = ITs[II];
                                //ss2 = "fit_N_" + N.ToString() + "_Pop_" + NPop.ToString() + "_intensity" + intensity.ToString() + "alpha_" + alpha.ToString() + "_c_" + c.ToString() + method + "_Acc_" + Ac.Text + "JR" + Jr.ToString() + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_struct_" + structure + "OF_" + OF.ToString() + ".txt";
                                //string ss1 = "JR_N_" + N.ToString() + "_Pop_" + NPop.ToString() + "_intensity" + intensity.ToString() + "alpha_" + alpha.ToString() + "_c_" + c.ToString() + method + "_Acc_" + Ac.Text + "JR" + Jr.ToString() + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_struct_" + structure + "OF_" + OF.ToString() + ".txt";
                                //string ss3 = "Rho_N_" + N.ToString() + "_Pop_" + NPop.ToString() + "_intensity" + intensity.ToString() + "alpha_" + alpha.ToString() + "_c_" + c.ToString() + method + "_Acc_" + Ac.Text + "JR" + Jr.ToString() + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_struct_" + structure + "OF_" + OF.ToString() + ".txt";
                                //string ss4 = "Alph_N_" + N.ToString() + "_Pop_" + NPop.ToString() + "_intensity" + intensity.ToString() + "alpha_" + alpha.ToString() + "_c_" + c.ToString() + method + "_Acc_" + Ac.Text + "JR" + Jr.ToString() + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_struct_" + structure + "OF_" + OF.ToString() + ".txt";
                                //string ss5 = "Fits_" + N.ToString() + "_Pop_" + NPop.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_struct_" + structure + "OF_" + OF.ToString() + ".txt";
                                //if (!File.Exists(ss5))
                                //{
                                //    for (int tt = 0; tt < Ntest; tt++)
                                //    {
                                //        int ik = 0;
                                //        IntenI = new double[NPop];
                                //        AlphaI = new double[NPop];
                                //        JrI = new double[NPop];
                                //        IntenAr = new ArrayList();
                                //        AlphaAR = new ArrayList();
                                //        JrAR = new ArrayList();
                                //        NF[] Nodes = new NF[NPop];
                                //        NF[] ONodes = new NF[NPop];
                                //        NF[] Forces = new NF[NPop];
                                //        NF[] a = new NF[NPop];
                                //        NF[] v = new NF[NPop];

                                //        intensity = 0.5;
                                //        alpha = 0.5;
                                //        Jr = 1.0;
                                //        double[] Magnet = new double[NPop];
                                //        double[] Mass = new double[NPop];
                                //        NF Best_indi = new NF();
                                //        double[] distances = new double[NPop];
                                //        Best_indi.F = -1E290;
                                //        for (int ii = 0; ii < NPop; ii++)
                                //        {
                                //            Forces[ii].x = new double[N];
                                //            a[ii].x = new double[N];
                                //            v[ii].x = new double[N];
                                //        }
                                //        for (int ii = 0; ii < NPop; ii++)
                                //        {
                                //            Nodes[ii].x = NFC.node_constructing(N, Rand);
                                //            //ONodes[ii].x = NFC.Opposition_Learning(Nodes[ii].x);
                                //            if (OF < 16 || OF == 22)
                                //                Nodes[ii].F = NFC.compute_fitness(Nodes[ii].x, OF, N);
                                //            else
                                //                Nodes[ii].F = NFC.compute_fitness(Nodes[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                //            //PervFitness[ii] = Nodes[ii].F;
                                //            //if (OF < 16 || OF == 22)
                                //            //    ONodes[ii].F = NFC.compute_fitness(ONodes[ii].x, OF, N);
                                //            //else
                                //            //    ONodes[ii].F = NFC.compute_fitness(ONodes[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                //            //if (Nodes[ii].F < ONodes[ii].F)
                                //            //{
                                //            //    Nodes[ii].x = (double[])ONodes[ii].x.Clone();
                                //            //    Nodes[ii].F = ONodes[ii].F;
                                //            //}
                                //        }
                                //        for (int it = 0; it < Iteration; it++)
                                //        {
                                //            IntenAr.Clear(); JrAR.Clear(); AlphaAR.Clear();
                                //            Fmax = -1E290; Fmin = 1000000000;
                                //            int min_index = 0; int max_index = 0;
                                //            for (int ii = 0; ii < NPop; ii++)
                                //            {
                                //                //JrI[ii] = NFC.Normal_Distribution(Rand, Jr, 0.5);
                                //                //RandGenNormal rg = new RandGenNormal(Jr, 0.01);
                                //                JrI[ii] = rg.NextDouble();
                                //                PervFitness[ii] = Nodes[ii].F;
                                //                if (OF < 16 || OF == 22)
                                //                    Nodes[ii].F = NFC.compute_fitness(Nodes[ii].x, OF, N);
                                //                else
                                //                    Nodes[ii].F = NFC.compute_fitness(Nodes[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                //                if (it > 0)
                                //                    if (Nodes[ii].F > PervFitness[ii])
                                //                    {
                                //                        IntenAr.Add(IntenI[ii]);
                                //                        AlphaAR.Add(AlphaI[ii]);

                                //                        ////////////Newly added
                                //                        JrAR.Add(JrI[ii]);
                                //                    }
                                //                if (JrI[ii] > Rand.NextDouble())
                                //                {

                                //                    ONodes[ii].x = NFC.generating_jumping(Nodes, NPop, N, ii);
                                //                    if (OF < 16 || OF == 22)
                                //                        ONodes[ii].F = NFC.compute_fitness(ONodes[ii].x, OF, N);
                                //                    else
                                //                        ONodes[ii].F = NFC.compute_fitness(ONodes[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                //                    if (Nodes[ii].F < ONodes[ii].F)
                                //                    {
                                //                        Nodes[ii].x = (double[])ONodes[ii].x.Clone();
                                //                        Nodes[ii].F = ONodes[ii].F;

                                //                    }
                                //                    Nodes[ii].x = NFC.Node_treatment(Nodes[ii].x, N, Rand);
                                //                }
                                //                if (Fmax < Nodes[ii].F)
                                //                {
                                //                    Fmax = Nodes[ii].F;
                                //                    max_index = ii;
                                //                }
                                //                if (Fmin > Nodes[ii].F)
                                //                {
                                //                    Fmin = Nodes[ii].F;
                                //                    min_index = ii;
                                //                }
                                //            }
                                //            ////////////////////////If the obl causes better values for the Indis the corresponding value is added.
                                //            //for (int jjs = 0; jjs < NPop; jjs++)
                                //            //{
                                //            //    if (Nodes[jjs].F < ONodes[jjs].F)
                                //            //        JrAR.Add(JrI[jjs]);
                                //            //}
                                //            if (Best_indi.F < Fmax)
                                //            {
                                //                Best_indi.F = Fmax;
                                //                Best_indi.x = (double[])Nodes[max_index].x.Clone();
                                //            }
                                //            double range = Fmax - Fmin;
                                //            if (range == 0)
                                //            {
                                //                Fmax = -1000000000; Fmin = 1000000000;
                                //                min_index = 0; max_index = 0;
                                //                for (int ii = 0; ii < NPop; ii++)
                                //                {
                                //                    Nodes[ii].x = NFC.node_constructing(N, Rand);
                                //                    ONodes[ii].x = NFC.Opposition_Learning(Nodes[ii].x);
                                //                }
                                //            }
                                //            for (int ii = 0; ii < NPop; ii++)
                                //            {
                                //                //IntenI[ii] = NFC.Normal_Distribution(Rand, intensity, 1);
                                //                //AlphaI[ii] = NFC.Normal_Distribution(Rand, alpha, 1);
                                //                RandGenNormal rg2 = new RandGenNormal(intensity, 0.1);
                                //                RandGenNormal rg3 = new RandGenNormal(alpha, 0.1);
                                //                IntenI[ii] = rg2.NextDouble();
                                //                AlphaI[ii] = rg3.NextDouble();
                                //                Magnet[ii] = (Nodes[ii].F - Fmin) / (range);
                                //                Mass[ii] = Magnet[ii] * IntenI[ii] + AlphaI[ii];
                                //            }
                                //            for (int ii = 0; ii < NPop; ii++)
                                //            {
                                //                Forces[ii].x = new double[N];
                                //            }
                                //            for (int ii = 0; ii < NPop; ii++)
                                //            {
                                //                int[] Neigbour = NFC.Structures(ii, NPop, structure);
                                //                for (int jj = 0; jj < Neigbour.Length; jj++)
                                //                {
                                //                    if (ii != Neigbour[jj])
                                //                    {
                                //                        if (D.Text == "1")
                                //                            distances[ii] = distance1(Nodes[ii].x, Nodes[Neigbour[jj]].x, N);
                                //                        //else
                                //                        //    distances[ii] = distance4(LHC[ii].LH, LHC[Neigbour[jj]].LH, N, P);
                                //                        if (distances[ii] > short_range)
                                //                        {
                                //                            for (int n = 0; n < N; n++)
                                //                                Forces[ii].x[n] = Forces[ii].x[n] + (Nodes[Neigbour[jj]].x[n] - Nodes[ii].x[n]) * Magnet[Neigbour[jj]] / distances[ii];

                                //                        }

                                //                    }
                                //                }
                                //                //for (int n = 0; n < N; n++)
                                //                //    for (int p = 0; p < P; p++)
                                //                //    {
                                //                //        Forces[ii].x[n, p] = Forces[ii].x[n, p] / (NPop - 1);
                                //                //        //Forces[ii].x[n, p] = Forces[ii].x[n, p] * rr.NextDouble();
                                //                //    }

                                //            }
                                //            if (acceleration == 0)
                                //            {
                                //                for (int ii = 0; ii < NPop; ii++)
                                //                {
                                //                    for (int n = 0; n < N; n++)
                                //                    {
                                //                        v[ii].x[n] = (Forces[ii].x[n] / Mass[ii]);
                                //                        Nodes[ii].x[n] = Nodes[ii].x[n] + v[ii].x[n];
                                //                    }

                                //                    Nodes[ii].x = NFC.Node_treatment(Nodes[ii].x, N, Rand);
                                //                }
                                //            }
                                //            else if (acceleration == 1)
                                //            {
                                //                for (int ii = 0; ii < NPop; ii++)
                                //                {
                                //                    for (int n = 0; n < N; n++)
                                //                    {
                                //                        a[ii].x[n] = (Forces[ii].x[n] / Mass[ii]);
                                //                        v[ii].x[n] = v[ii].x[n] * Rand.NextDouble() + a[ii].x[n];
                                //                        Nodes[ii].x[n] = Nodes[ii].x[n] + v[ii].x[n];
                                //                    }

                                //                    Nodes[ii].x = NFC.Node_treatment(Nodes[ii].x, N, Rand);
                                //                }
                                //            }
                                //            //////update///////////////////////////////////////
                                //            if (JrAR.Count > 0)
                                //            {
                                //                Jr = (1 - c) * Jr + c * NFC.mean_simple(JrAR);
                                //                //Jr = Rand.NextDouble() * 0.3 + 0.3;
                                //            }
                                //            //if (it > 0)
                                //            //{
                                //            if (AlphaAR.Count > 0)
                                //            {
                                //                alpha = (1 - c) * alpha + c * NFC.mean_Lehmer(AlphaAR);
                                //                //alpha = Rand.NextDouble() * 100;
                                //                //intensity = Rand.NextDouble() * 100;
                                //                intensity = (1 - c) * intensity + c * NFC.mean_Lehmer(IntenAr);
                                //            }
                                //            //}
                                //            ///////////////////////////////////////////////////
                                //            //if ((it % 100) == 1)
                                //            //{
                                //            Js[tt, ik] = Jr;
                                //            Alphas[tt, ik] = alpha;
                                //            Rhos[tt, ik] = intensity;
                                //            Fs[tt, ik++] = Best_indi.F;
                                //            //}
                                //        }
                                //        costs[tt] = Best_indi.F;
                                //        //Nodes_test[tt].x = (double[])Best_indi.x.Clone();
                                //    }
                                //    sw = new StreamWriter(ss2);
                                //    for (int tt = 0; tt < Ntest; tt++)
                                //    {
                                //        sw.Write(costs[tt].ToString() + " ");
                                //    }
                                //    sw.Close();
                                //    sw = new StreamWriter(ss1);
                                //    for (int tt = 0; tt < Ntest; tt++)
                                //    {
                                //        for (int it = 0; it < Iteration; it++)
                                //        {
                                //            sw.Write(Js[tt, it].ToString() + " ");
                                //        }
                                //        sw.WriteLine();
                                //    }
                                //    sw.Close();
                                //    sw = new StreamWriter(ss3);
                                //    for (int tt = 0; tt < Ntest; tt++)
                                //    {
                                //        for (int it = 0; it < Iteration; it++)
                                //        {
                                //            sw.Write(Rhos[tt, it].ToString() + " ");
                                //        }
                                //        sw.WriteLine();
                                //    }
                                //    sw.Close();
                                //    sw = new StreamWriter(ss4);
                                //    for (int tt = 0; tt < Ntest; tt++)
                                //    {
                                //        for (int it = 0; it < Iteration; it++)
                                //        {
                                //            sw.Write(Alphas[tt, it].ToString() + " ");
                                //        }
                                //        sw.WriteLine();
                                //    }
                                //    sw.Close();
                                //    sw = new StreamWriter(ss5);
                                //    for (int tt = 0; tt < Ntest; tt++)
                                //    {
                                //        for (int it = 0; it < Iteration; it++)
                                //        {
                                //            sw.Write(Fs[tt, it].ToString() + " ");
                                //        }
                                //        sw.WriteLine();
                                //    }
                                //    sw.Close();
                                //}
                            }
                            #endregion
                            #region MOA with the parameters set by iterated-Frace
                            if (method == "FMOAFR")
                            {
                                double Jr = Convert.ToDouble(Jumping_rate.Text);
                                double short_range = 0.01;
                                int acceleration = Convert.ToInt32(Ac.Text);
                                int Distance = Convert.ToInt32(D.Text);
                                double al1 = Convert.ToDouble(Al.Text);
                                double al2 = Convert.ToDouble(Al2.Text);
                                double al_step = Convert.ToDouble(Al_step.Text);
                                double in1 = Convert.ToDouble(Inten.Text);
                                double in2 = Convert.ToDouble(Inten2.Text);
                                double in_step = Convert.ToDouble(Inten_step.Text);
                                //double[] Intinsities = { 0.001, 0.01, 0.2, 0.6, 1 };
                                //double[] Alphas = { 0.001, 0.01, 0.2, 0.6, 1 };
                                //double[] Jrs = { 0.3, 0.45, 0.6 };
                                //double[] Intinsities = { 0.2, 0.001, 0.01, 0.2, 0.2, 0.01 };
                                //double[] Alphas = { 1, 0.6, 0.2, 0.001, 0.001, 0.2 };
                                //double[] Jrs = { 0.3, 0.45, 0.6, 0.6, 0.6, 0.6 };
                                //int ITT = 0;
                                //double[] Intinsities = { 0.6, 1.0, 1.0, 0.2, 0.6, 1.0, 0.2, 1.0, 1.0, 1.0, 1.0, 0.6, 1.0, 0, 1.0, 0.2, 0.2, 0.0010, 0.0010, 0.00100, 0.0010, 1.0 };
                                //double[] Alphas = { 0.0200000000000000, 0.400000000000000, 0.400000000000000, 0.200000000000000, 1, 0.800000000000000, 0.200000000000000, 0.200000000000000, 0.6, 0.4, 0.4, 1, 0.2, 0, 0.6, 1, 0.6, 0.4, 0.0010, 0.01000, 0.2000, 1 };
                                //double[] Jrs = { 0.6, 0.3, 0.45, 0.3, 0.3, 0.3, 0.45, 0.3, 0.3, 0.3, 0.3, 0.3, 0.45, 0, 0.3000, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.3 };
                                //int[] POP = { 5, 10, 25, 50, 100, 200, 250, 500 };
                                //int[] ITs = { 10000, 5000, 2000, 1000, 500, 250, 200, 100 };
                                //for (int al_index = 0; al_index < Alphas.Length; al_index++)
                                //    for (int int_index = 0; int_index < Intinsities.Length; int_index++)
                                //        for (int II = 0; II < 3; II++)
                                //        {
                                double intensity = 93.01561;
                                double alpha = 86.6102;
                                Jr = 0.7958435;
                                //double intensity = Intinsities[OF - 22];
                                //double alpha = Alphas[OF - 22];
                                //Jr = Jrs[OF - 22];
                                //NPop = POP[II];
                                //Iteration = ITs[II];
                                //ss = "Result_N_" + N.ToString() + "_Pop_" + NPop.ToString() + "_intensity" + intensity.ToString() + "_alpha_" + alpha.ToString() + "_Dis_" + method + Distance.ToString() + "_Acc_" + Ac.Text + "JR" + Jr.ToString() + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_struct_" + structure + "OF_" + OF.ToString() + ".txt";
                                //ss2 = "Results\\fit_N_" + N.ToString() + "_Pop_" + NPop.ToString() + "_intensity" + intensity.ToString() + "alpha_" + alpha.ToString() + "_Dis_" + method + Distance.ToString() + "_Acc_" + Ac.Text + "JR" + Jr.ToString() + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_struct_" + structure + "OF_" + OF.ToString() + ".txt";
                                if (!File.Exists(result))
                                {
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        NF[] Nodes = new NF[NPop];
                                        NF[] ONodes = new NF[NPop];
                                        NF[] Forces = new NF[NPop];
                                        NF[] a = new NF[NPop];
                                        NF[] v = new NF[NPop];
                                        double[] Magnet = new double[NPop];
                                        double[] Mass = new double[NPop];
                                        NF Best_indi = new NF();
                                        double[] distances = new double[NPop];
                                        Best_indi.F = -1E290;
                                        for (int ii = 0; ii < NPop; ii++)
                                        {
                                            Forces[ii].x = new double[N];
                                            a[ii].x = new double[N];
                                            v[ii].x = new double[N];
                                        }
                                        for (int ii = 0; ii < NPop; ii++)
                                        {
                                            Nodes[ii].x = NFC.node_constructing(N, Rand);
                                            ONodes[ii].x = NFC.Opposition_Learning(Nodes[ii].x);
                                        }
                                        for (int it = 0; it < Iteration; it++)
                                        {
                                            Fmax = -1E290; Fmin = 1000000000;
                                            int min_index = 0; int max_index = 0;
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                if (OF < 16 || OF == 22)
                                                    Nodes[ii].F = NFC.compute_fitness(Nodes[ii].x, OF, N);
                                                else
                                                    Nodes[ii].F = NFC.compute_fitness(Nodes[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);

                                                if (Jr > Rand.NextDouble())
                                                {

                                                    ONodes[ii].x = NFC.generating_jumping(Nodes, NPop, N, ii);
                                                    if (OF < 16 || OF == 22)
                                                        ONodes[ii].F = NFC.compute_fitness(ONodes[ii].x, OF, N);
                                                    else
                                                        ONodes[ii].F = NFC.compute_fitness(ONodes[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                    if (Nodes[ii].F < ONodes[ii].F)
                                                    {
                                                        Nodes[ii].x = (double[])ONodes[ii].x.Clone();
                                                        Nodes[ii].F = ONodes[ii].F;
                                                    }
                                                    Nodes[ii].x = NFC.Node_treatment(Nodes[ii].x, N, Rand);
                                                }

                                                if (Fmax < Nodes[ii].F)
                                                {
                                                    Fmax = Nodes[ii].F;
                                                    max_index = ii;
                                                }
                                                if (Fmin > Nodes[ii].F)
                                                {
                                                    Fmin = Nodes[ii].F;
                                                    min_index = ii;
                                                }
                                            }

                                            if (Best_indi.F < Fmax)
                                            {
                                                Best_indi.F = Fmax;
                                                Best_indi.x = (double[])Nodes[max_index].x.Clone();
                                            }
                                            double range = Fmax - Fmin;
                                            if (range == 0)
                                            {
                                                Fmax = -1000000000; Fmin = 1000000000;
                                                min_index = 0; max_index = 0;
                                                for (int ii = 0; ii < NPop; ii++)
                                                {
                                                    Nodes[ii].x = NFC.node_constructing(N, Rand);
                                                    ONodes[ii].x = NFC.Opposition_Learning(Nodes[ii].x);
                                                }
                                            }
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                Magnet[ii] = (Nodes[ii].F - Fmin) / (range);
                                                Mass[ii] = Magnet[ii] * intensity + alpha;
                                            }
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                Forces[ii].x = new double[N];
                                            }
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                int[] Neigbour = NFC.Structures(ii, NPop, structure);
                                                for (int jj = 0; jj < Neigbour.Length; jj++)
                                                {
                                                    if (ii != Neigbour[jj])
                                                    {
                                                        if (D.Text == "1")
                                                            distances[ii] = distance1(Nodes[ii].x, Nodes[Neigbour[jj]].x, N);
                                                        //else
                                                        //    distances[ii] = distance4(LHC[ii].LH, LHC[Neigbour[jj]].LH, N, P);
                                                        if (distances[ii] > short_range)
                                                        {
                                                            for (int n = 0; n < N; n++)
                                                                Forces[ii].x[n] = Forces[ii].x[n] + (Nodes[Neigbour[jj]].x[n] - Nodes[ii].x[n]) * Magnet[Neigbour[jj]] / distances[ii];

                                                        }

                                                    }
                                                }
                                                //for (int n = 0; n < N; n++)
                                                //    for (int p = 0; p < P; p++)
                                                //    {
                                                //        Forces[ii].x[n, p] = Forces[ii].x[n, p] / (NPop - 1);
                                                //        //Forces[ii].x[n, p] = Forces[ii].x[n, p] * rr.NextDouble();
                                                //    }

                                            }
                                            if (acceleration == 0)
                                            {
                                                for (int ii = 0; ii < NPop; ii++)
                                                {
                                                    for (int n = 0; n < N; n++)
                                                    {
                                                        v[ii].x[n] = (Forces[ii].x[n] / Mass[ii]);
                                                        Nodes[ii].x[n] = Nodes[ii].x[n] + v[ii].x[n];
                                                    }

                                                    Nodes[ii].x = NFC.Node_treatment(Nodes[ii].x, N, Rand);
                                                }
                                            }
                                            else if (acceleration == 1)
                                            {
                                                for (int ii = 0; ii < NPop; ii++)
                                                {
                                                    for (int n = 0; n < N; n++)
                                                    {
                                                        a[ii].x[n] = (Forces[ii].x[n] / Mass[ii]);
                                                        v[ii].x[n] = v[ii].x[n] * Rand.NextDouble() + a[ii].x[n];
                                                        Nodes[ii].x[n] = Nodes[ii].x[n] + v[ii].x[n];
                                                    }

                                                    Nodes[ii].x = NFC.Node_treatment(Nodes[ii].x, N, Rand);
                                                }
                                            }


                                        }
                                        costs[tt] = Best_indi.F;
                                        Nodes_test[tt].x = (double[])Best_indi.x.Clone();
                                    }
                                    sw = new StreamWriter(result);
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        sw.Write(costs[tt].ToString() + " ");
                                    }
                                    sw.Close();
                                }

                            }
                            #endregion
                            #region Analaysis 2
                            if (method == "FMOAAN2")
                            {
                                double Jr = Convert.ToDouble(Jumping_rate.Text);
                                double short_range = 0.01;
                                int acceleration = Convert.ToInt32(Ac.Text);
                                int Distance = Convert.ToInt32(D.Text);
                                double al1 = Convert.ToDouble(Al.Text);
                                double al2 = Convert.ToDouble(Al2.Text);
                                double al_step = Convert.ToDouble(Al_step.Text);
                                double in1 = Convert.ToDouble(Inten.Text);
                                double in2 = Convert.ToDouble(Inten2.Text);
                                double in_step = Convert.ToDouble(Inten_step.Text);
                                double[] Intinsities = { 0.001, 0.01, 0.2, 0.6, 1, 2, 5, 10, 30, 50, 100 };
                                double[] Alphas = { 0.001, 0.01, 0.2, 0.6, 1, 2, 5, 10, 30, 50, 100 };
                                double[] Jrs = { 0.3, 0.45, 0.6, 0.8, 1 };
                                //double[] Intinsities = { 0.2, 0.001, 0.01, 0.2, 0.2, 0.01 };
                                //double[] Alphas = { 1, 0.6, 0.2, 0.001, 0.001, 0.2 };
                                //double[] Jrs = { 0.3, 0.45, 0.6, 0.6, 0.6, 0.6 };
                                //int ITT = 0;
                                double[,] fits = new double[Ntest, Iteration];
                                //double[] Intinsities = { 0.6, 1.0, 1.0, 0.2, 0.6, 1.0, 0.2, 1.0, 1.0, 1.0, 1.0, 0.6, 1.0, 0, 1.0, 0.2, 0.2, 0.0010, 0.0010, 0.00100, 0.0010, 1.0 };
                                //double[] Alphas = { 0.0200000000000000, 0.400000000000000, 0.400000000000000, 0.200000000000000, 1, 0.800000000000000, 0.200000000000000, 0.200000000000000, 0.6, 0.4, 0.4, 1, 0.2, 0, 0.6, 1, 0.6, 0.4, 0.0010, 0.01000, 0.2000, 1 };
                                //double[] Jrs = { 0.6, 0.3, 0.45, 0.3, 0.3, 0.3, 0.45, 0.3, 0.3, 0.3, 0.3, 0.3, 0.45, 0, 0.3000, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.3 };
                                //int[] POP = { 5, 10, 25, 50, 100, 200, 250, 500 };
                                //int[] ITs = { 10000, 5000, 2000, 1000, 500, 250, 200, 100 };
                                //for (int al_index = 0; al_index < Alphas.Length; al_index++)
                                //    for (int int_index = 0; int_index < Intinsities.Length; int_index++)
                                //        for (int II = 0; II < 3; II++)
                                //        {
                                //double intensity = 93.01561;
                                //double alpha = 86.6102;
                                //Jr = 0.7958435;
                                NF[, , ,] TempXs = new NF[Intinsities.Length, Alphas.Length, Jrs.Length, NPop];
                                NF[, , ,] TempOXs = new NF[Intinsities.Length, Alphas.Length, Jrs.Length, NPop];
                                double[, ,] par_indexes = new double[Ntest, Iteration, 3];
                                int par_in1 = -1;
                                int par_in2 = -1;
                                int par_in3 = -1;
                                NF[, ,] Best_indis = new NF[Intinsities.Length, Alphas.Length, Jrs.Length];
                                double fit_par = -1;
                                //NPop = POP[II];
                                //Iteration = ITs[II];
                                //ss = "Result_N_" + N.ToString() + "_Pop_" + NPop.ToString() + "_intensity" + intensity.ToString() + "_alpha_" + alpha.ToString() + "_Dis_" + method + Distance.ToString() + "_Acc_" + Ac.Text + "JR" + Jr.ToString() + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_struct_" + structure + "OF_" + OF.ToString() + ".txt";
                                //ss2 = "fit_N_" + N.ToString() + "_Pop_" + NPop.ToString() + "_Dis_" + method + Distance.ToString() + "_Acc_" + Ac.Text + "JR" + Jr.ToString() + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_struct_" + structure + "OF_" + OF.ToString() + ".txt";
                                //string ss3 = "fits_N_" + N.ToString() + "_Pop_" + NPop.ToString() + "_Dis_" + method + Distance.ToString() + "_Acc_" + Ac.Text + "JR" + Jr.ToString() + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_struct_" + structure + "OF_" + OF.ToString() + ".txt";
                                //string par_val = "par_val_N_" + N.ToString() + "_Pop_" + NPop.ToString() + "_Dis_" + method + Distance.ToString() + "_Acc_" + Ac.Text + "JR" + Jr.ToString() + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_struct_" + structure + "OF_" + OF.ToString() + ".txt";
                                //string par_val2 = "par_val2_N_" + N.ToString() + "_Pop_" + NPop.ToString() + "_Dis_" + method + Distance.ToString() + "_Acc_" + Ac.Text + "JR" + Jr.ToString() + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_struct_" + structure + "OF_" + OF.ToString() + ".txt";
                                //string par_val3 = "par_val3_N_" + N.ToString() + "_Pop_" + NPop.ToString() + "_Dis_" + method + Distance.ToString() + "_Acc_" + Ac.Text + "JR" + Jr.ToString() + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_struct_" + structure + "OF_" + OF.ToString() + ".txt";
                                if (!File.Exists(result))
                                {
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        NF[] Nodes = new NF[NPop];
                                        NF[] ONodes = new NF[NPop];
                                        NF[] Forces = new NF[NPop];
                                        NF[] a = new NF[NPop];
                                        NF[] v = new NF[NPop];
                                        double[] Magnet = new double[NPop];
                                        double[] Mass = new double[NPop];
                                        NF Best_indi = new NF();
                                        double[] distances = new double[NPop];
                                        Best_indi.F = -1E290;
                                        for (int ii = 0; ii < NPop; ii++)
                                        {
                                            Forces[ii].x = new double[N];
                                            a[ii].x = new double[N];
                                            v[ii].x = new double[N];
                                        }
                                        for (int ii = 0; ii < NPop; ii++)
                                        {
                                            Nodes[ii].x = NFC.node_constructing(N, Rand);
                                            ONodes[ii].x = NFC.Opposition_Learning(Nodes[ii].x);
                                        }
                                        fit_par = -1E300;
                                        par_in1 = -1;
                                        par_in2 = -1;
                                        par_in3 = -1;
                                        Best_indi.F = -1E300;
                                        for (int it = 0; it < Iteration; it++)
                                        {

                                            for (int ht = 0; ht < Intinsities.Length; ht++)
                                                for (int ht2 = 0; ht2 < Alphas.Length; ht2++)
                                                    for (int ht3 = 0; ht3 < Jrs.Length; ht3++)
                                                    {
                                                        Fmax = -1E300; Fmin = 1000000000;
                                                        int min_index = 0; int max_index = 0;
                                                        if (it == 0)
                                                        {
                                                            //Best_indis[ht, ht2, ht3].x = (double[])Best_indi.x.Clone();
                                                        }
                                                        else
                                                        {
                                                            Best_indi.x = (double[])Best_indis[ht, ht2, ht3].x.Clone();
                                                            Best_indi.F = Best_indis[ht, ht2, ht3].F;
                                                        }
                                                        for (int ii = 0; ii < NPop; ii++)
                                                        {

                                                            if (it == 0)
                                                            {
                                                                TempOXs[ht, ht2, ht3, ii].x = (double[])ONodes[ii].x.Clone();
                                                                TempXs[ht, ht2, ht3, ii].x = (double[])Nodes[ii].x.Clone();
                                                            }
                                                            else
                                                            {
                                                                ONodes[ii].x = (double[])TempOXs[ht, ht2, ht3, ii].x.Clone();
                                                                Nodes[ii].x = (double[])TempXs[ht, ht2, ht3, ii].x.Clone();
                                                            }
                                                            if (OF < 16 || OF == 22)
                                                                Nodes[ii].F = NFC.compute_fitness(Nodes[ii].x, OF, N);
                                                            else
                                                                Nodes[ii].F = NFC.compute_fitness(Nodes[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);

                                                            if (Jrs[ht3] > Rand.NextDouble())
                                                            {

                                                                ONodes[ii].x = NFC.generating_jumping(Nodes, NPop, N, ii);
                                                                if (OF < 16 || OF == 22)
                                                                    ONodes[ii].F = NFC.compute_fitness(ONodes[ii].x, OF, N);
                                                                else
                                                                    ONodes[ii].F = NFC.compute_fitness(ONodes[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                                if (Nodes[ii].F < ONodes[ii].F)
                                                                {
                                                                    Nodes[ii].x = (double[])ONodes[ii].x.Clone();
                                                                    Nodes[ii].F = ONodes[ii].F;

                                                                }
                                                                Nodes[ii].x = NFC.Node_treatment(Nodes[ii].x, N, Rand);
                                                            }

                                                            if (Fmax < Nodes[ii].F)
                                                            {
                                                                Fmax = Nodes[ii].F;
                                                                max_index = ii;
                                                            }
                                                            if (Fmin > Nodes[ii].F)
                                                            {
                                                                Fmin = Nodes[ii].F;
                                                                min_index = ii;
                                                            }
                                                        }

                                                        if (Best_indi.F < Fmax)
                                                        {
                                                            Best_indi.F = Fmax;
                                                            Best_indi.x = (double[])Nodes[max_index].x.Clone();
                                                        }
                                                        double range = Fmax - Fmin;
                                                        if (range == 0)
                                                        {
                                                            Fmax = -1000000000; Fmin = 1000000000;
                                                            min_index = 0; max_index = 0;
                                                            for (int ii = 0; ii < NPop; ii++)
                                                            {
                                                                Nodes[ii].x = NFC.node_constructing(N, Rand);
                                                                ONodes[ii].x = NFC.Opposition_Learning(Nodes[ii].x);
                                                            }
                                                        }
                                                        for (int ii = 0; ii < NPop; ii++)
                                                        {
                                                            Magnet[ii] = (Nodes[ii].F - Fmin) / (range);
                                                            Mass[ii] = Magnet[ii] * Intinsities[ht] + Alphas[ht2];
                                                        }
                                                        for (int ii = 0; ii < NPop; ii++)
                                                        {
                                                            Forces[ii].x = new double[N];
                                                        }
                                                        for (int ii = 0; ii < NPop; ii++)
                                                        {
                                                            int[] Neigbour = NFC.Structures(ii, NPop, structure);
                                                            for (int jj = 0; jj < Neigbour.Length; jj++)
                                                            {
                                                                if (ii != Neigbour[jj])
                                                                {
                                                                    if (D.Text == "1")
                                                                        distances[ii] = distance1(Nodes[ii].x, Nodes[Neigbour[jj]].x, N);
                                                                    //else
                                                                    //    distances[ii] = distance4(LHC[ii].LH, LHC[Neigbour[jj]].LH, N, P);
                                                                    if (distances[ii] > short_range)
                                                                    {
                                                                        for (int n = 0; n < N; n++)
                                                                            Forces[ii].x[n] = Forces[ii].x[n] + (Nodes[Neigbour[jj]].x[n] - Nodes[ii].x[n]) * Magnet[Neigbour[jj]] / distances[ii];

                                                                    }

                                                                }
                                                            }
                                                            //for (int n = 0; n < N; n++)
                                                            //    for (int p = 0; p < P; p++)
                                                            //    {
                                                            //        Forces[ii].x[n, p] = Forces[ii].x[n, p] / (NPop - 1);
                                                            //        //Forces[ii].x[n, p] = Forces[ii].x[n, p] * rr.NextDouble();
                                                            //    }

                                                        }
                                                        if (acceleration == 0)
                                                        {
                                                            for (int ii = 0; ii < NPop; ii++)
                                                            {
                                                                for (int n = 0; n < N; n++)
                                                                {
                                                                    v[ii].x[n] = (Forces[ii].x[n] / Mass[ii]);
                                                                    Nodes[ii].x[n] = Nodes[ii].x[n] + v[ii].x[n];
                                                                }

                                                                Nodes[ii].x = NFC.Node_treatment(Nodes[ii].x, N, Rand);
                                                            }
                                                        }
                                                        else if (acceleration == 1)
                                                        {
                                                            for (int ii = 0; ii < NPop; ii++)
                                                            {
                                                                for (int n = 0; n < N; n++)
                                                                {
                                                                    a[ii].x[n] = (Forces[ii].x[n] / Mass[ii]);
                                                                    v[ii].x[n] = v[ii].x[n] * Rand.NextDouble() + a[ii].x[n];
                                                                    Nodes[ii].x[n] = Nodes[ii].x[n] + v[ii].x[n];
                                                                }

                                                                Nodes[ii].x = NFC.Node_treatment(Nodes[ii].x, N, Rand);
                                                            }
                                                        }
                                                        if (Best_indi.F > fit_par)
                                                        {
                                                            par_in1 = ht;
                                                            par_in2 = ht2;
                                                            par_in3 = ht3;
                                                            fit_par = Best_indi.F;
                                                        }
                                                        for (int ii = 0; ii < NPop; ii++)
                                                        {
                                                            TempOXs[ht, ht2, ht3, ii].x = (double[])ONodes[ii].x.Clone();
                                                            TempXs[ht, ht2, ht3, ii].x = (double[])Nodes[ii].x.Clone();
                                                            TempXs[ht, ht2, ht3, ii].F = Nodes[ii].F;
                                                            TempOXs[ht, ht2, ht3, ii].F = ONodes[ii].F;
                                                        }
                                                        Best_indis[ht, ht2, ht3].x = (double[])Best_indi.x.Clone();
                                                        Best_indis[ht, ht2, ht3].F = Best_indi.F;
                                                    }
                                            par_indexes[tt, it, 0] = Intinsities[par_in1];
                                            par_indexes[tt, it, 1] = Alphas[par_in2];
                                            par_indexes[tt, it, 2] = Jrs[par_in3];
                                            fits[tt, it] = fit_par;
                                        }
                                    }
                                    sw = new StreamWriter(result);
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        for (int ii = 0; ii < Iteration; ii++)
                                        {
                                            sw.Write(fits[tt, ii].ToString() + " ");
                                        }
                                        sw.WriteLine();
                                    }
                                    sw.Close();
                                    //sw = new StreamWriter(par_val);
                                    //for (int tt = 0; tt < Ntest; tt++)
                                    //{
                                    //    for (int ii = 0; ii < Iteration; ii++)
                                    //    {
                                    //        sw.Write(par_indexes[tt, ii, 0].ToString() + " ");
                                    //    }
                                    //    sw.WriteLine();
                                    //}
                                    //sw.Close();
                                    //sw = new StreamWriter(par_val2);
                                    //for (int tt = 0; tt < Ntest; tt++)
                                    //{
                                    //    for (int ii = 0; ii < Iteration; ii++)
                                    //    {
                                    //        sw.Write(par_indexes[tt, ii, 1].ToString() + " ");
                                    //    }
                                    //    sw.WriteLine();
                                    //}
                                    //sw.Close();
                                    //sw = new StreamWriter(par_val3);
                                    //for (int tt = 0; tt < Ntest; tt++)
                                    //{
                                    //    for (int ii = 0; ii < Iteration; ii++)
                                    //    {
                                    //        sw.Write(par_indexes[tt, ii, 2].ToString() + " ");
                                    //    }
                                    //    sw.WriteLine();
                                    //}
                                    //sw.Close();
                                }
                            }
                            #endregion
                            #region CMA-ES
                            if (method == "CMAES")
                            {
                                double sigma = 0.3;// is proposed by the wiki
                                //xmean = rand(N,1);    % objective variables initial point
                                //sigma = 0.3;          % coordinate wise standard deviation (step size)
                                int lambda = 4 + Convert.ToInt32(Math.Floor(3 * Math.Log(N)));
                                int mu = lambda / 2;
                                double mueff = -1;
                                int counteval = 0;
                                Random rr = new Random();
                                double[,] xmean = new double[N, 1];
                                for (int ii = 0; ii < N; ii++)
                                {
                                    xmean[ii, 0] = rr.NextDouble();
                                }
                                // Strategy parameter setting: Selection  
                                double[,] weights = new double[mu, 1];
                                double sum = 0;
                                double sum2 = 0;
                                for (int ii = 0; ii < mu; ii++)
                                {
                                    weights[ii, 0] = Math.Log(mu + 1 / 2) - Math.Log(ii + 1);
                                    sum += weights[ii, 0];
                                    sum2 += Math.Pow(weights[ii, 0], 2);
                                }
                                for (int ii = 0; ii < mu; ii++)
                                {
                                    weights[ii, 0] = weights[ii, 0] / sum;
                                }
                                mueff = Math.Pow(sum, 2) / sum2;
                                // Strategy parameter setting: Adaptation
                                double cc = (4 + mueff / N) / (N + 4 + 2 * mueff / N); // time constant for cumulation for C
                                double cs = (mueff + 2) / (N + mueff + 5); // % t-const for cumulation for sigma control
                                double c1 = 2 / (Math.Pow((N + 1.3), 2) + mueff);//    % learning rate for rank-one update of C
                                double cmu = Math.Min(1 - c1, 2 * (mueff - 2 + 1 / mueff) / (Math.Pow((N + 2), 2) + mueff));//  % and for rank-mu update
                                double damps = 1 + 2 * Math.Max(0, Math.Sqrt((mueff - 1) / (N + 1)) - 1) + cs; //% damping for sigma % usually close to 1
                                //% Initialize dynamic (internal) strategy parameters and constants
                                double[] pc = new double[N];
                                double[] ps = new double[N];//   % evolution paths for C and sigma
                                double[,] B = function_evaluation2013.Eye(N);//    % B defines the coordinate system
                                double[,] D = new double[N, 1];
                                for (int ii = 0; ii < N; ii++)/// % diagonal D defines the scaling
                                {
                                    D[ii, 0] = 1;
                                }
                                // C = B * diag(D.^2) * B';           // % covariance matrix C
                                double[,] C = new double[N, N];
                                C = function_evaluation2013.Eye(N);
                                function_evaluation2013 fnewclass = new function_evaluation2013();
                                double[,] invsqrtC = fnewclass.MultiplyMatrix(fnewclass.MultiplyMatrix(B, fnewclass.diags(D)), fnewclass.replace2d(B));  //  % C^-1/2 
                                double eigeneval = 0;    // % track update of B and D
                                double chiN = Math.Pow(N, 0.5) * (1 - 1 / (4 * N) + 1 / (21 * N ^ 2));//  % expectation of ||N(0,I)|| == norm(randn(N,1))
                                //  % -------------------- Generation Loop --------------------------------
                                //ss2 = "Results\\fit_N_" + N.ToString() + "_Pop_" + NPop.ToString() + "_sigma_" + sigma.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "OF_" + OF.ToString() + ".txt";
                                if (!File.Exists(result))
                                {
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {

                                        NF[] Nodes = new NF[NPop];
                                        NF Best_nodes = new NF();
                                        Best_nodes.F = -1E250;
                                        int couneval = 0;
                                        int countsum = 50000;
                                        for (int ii = 0; ii < lambda; ii++)
                                        {
                                            Nodes[ii].x = NFC.node_constructing(N, Rand);

                                            if (OF < 16 || OF == 22)
                                                Nodes[ii].F = NFC.compute_fitness(Nodes[ii].x, OF, N);
                                            else
                                                Nodes[ii].F = NFC.compute_fitness(Nodes[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);

                                        }
                                        ///////////////////////////////////////////////////////////////////////////////////////
                                        while (couneval < countsum)
                                        {

                                            for (int ii = 0; ii < lambda; ii++)
                                            {

                                                // % Generate and evaluate lambda offspring
                                                Nodes[ii].x = fnewclass.create_new_gen_cma_es(N, xmean, sigma, B, D);  //% m + sig * Normal(0,C) 

                                                if (OF < 16 || OF == 22)
                                                    Nodes[ii].F = NFC.compute_fitness(Nodes[ii].x, OF, N);
                                                else
                                                    Nodes[ii].F = NFC.compute_fitness(Nodes[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                counteval = counteval + 1;
                                            }
                                            double[] arfitness = new double[NPop];
                                            int[] arindex = new int[NPop];
                                            Numerical_functions.Sort(Nodes, ref arindex, ref arfitness);
                                            //[arfitness, arindex] = sort(arfitness); //% minimization
                                            double[,] xold = (double[,])xmean.Clone();
                                            ////////////Compute x-related weighted matrix
                                            double[,] popselected = new double[N, mu];
                                            for (int ii = 0; ii < mu; ii++)
                                                for (int jjj = 0; jjj < N; jjj++)
                                                {
                                                    popselected[jjj, ii] = Nodes[arindex[ii]].x[jjj];
                                                }
                                            xmean = fnewclass.MultiplyMatrix(popselected, weights);


                                            //% Cumulation: Update evolution paths
                                            double hsig = fnewclass.Update_evolution_paths(cs, ref ps, mueff, invsqrtC, xmean, xold, sigma, couneval, lambda, chiN, N, cc, ref pc);
                                            // % Adapt covariance matrix C
                                            double[,] artmp = new double[mu, N];
                                            for (int ii = 0; ii < mu; ii++)
                                                for (int jj = 0; jj < N; jj++)
                                                {
                                                    artmp[ii, jj] = (1 / sigma) * popselected[jj, ii] - xold[jj, 0];
                                                }

                                            // cmu * artmp[ii,jj] *fnewclass.diags(weights) * artmp'
                                            double[,] wei = fnewclass.diags(weights);
                                            double[,] artmpsum = new double[artmp.GetLength(0), artmp.GetLength(0)];
                                            double[,] artmpapp = fnewclass.replace2d(artmp);
                                            double[,] tem = fnewclass.MultiplyMatrix(artmp, wei);
                                            double[,] tem2 = fnewclass.MultiplyMatrix(tem, artmpapp);
                                            for (int ii = 0; ii < artmpapp.GetLength(0); ii++)
                                                for (int jj = 0; jj < artmpapp.GetLength(1); jj++)
                                                    artmpsum[ii, jj] = cmu * tem2[ii, jj];

                                            for (int ii = 0; ii < C.GetLength(0); ii++)
                                                for (int jj = 0; jj < C.GetLength(1); jj++)
                                                {
                                                    C[ii, jj] = (1 - c1 - cmu) * C[ii, jj] + c1 * (fnewclass.replace(pc) + (1 - hsig) * cc * (2 - cc) * C[ii, jj]) + artmpsum[ii, jj];
                                                }

                                            //// % Adapt step size sigma
                                            ////////////////norm ps
                                            double norm = 0;
                                            for (int ii = 0; ii < ps.Length; ii++)
                                            {
                                                norm += Math.Pow(ps[ii], 2);
                                            }
                                            norm = Math.Sqrt(norm);
                                            sigma = sigma * Math.Exp((cs / damps) * (norm / chiN - 1));
                                            //% Decomposition of C into B*diag(D.^2)*B' (diagonalization)
                                            if (counteval - eigeneval > lambda / (c1 + cmu) / N / 10)//  % to achieve O(N^2)
                                            {
                                                //eigeneval = counteval;
                                                //C = triu(C) + triu(C,1)';// % enforce symmetry
                                                //[B,D] = eig(C);          // % eigen decomposition, B==normalized eigenvectors
                                                //D = sqrt(diag(D));        //% D is a vector of standard deviations now
                                                //invsqrtC = B * diag(D.^-1) * B';
                                            }

                                            //% Break, if fitness is good enough or condition exceeds 1e14, better termination methods are advisable 
                                        }

                                        //      % Sort by fitness and compute weighted mean into xmean





                                        costs[tt] = Best_nodes.F;
                                        Range.Refresh();
                                    }

                                    sw = new StreamWriter(result);
                                    for (int tt = 0; tt < Ntest; tt++)
                                    {
                                        sw.Write(costs[tt].ToString() + " ");
                                    }
                                    sw.Close();
                                }
                            }
                            #endregion
                        }
                    }
                }
            }
                #endregion
                #region RLD
            else if (Criterion.Text == "RLD")
            {
                for (int N = ND1; N <= ND2; N++)
                    for (int OF = OF1; OF <= OF2; OF += OF_s)
                    {
                        if (OF != 14)
                        {
                            Random rrr = new Random(N);
                            int M = 20;
                            int[] Perm = NFC.sq_rand_gen2(N, N, rrr);
                            double[,] Orthogonal_Matrix;
                            double[] o = new double[N];
                            Random RO = new Random(100);
                            o = NFC.node_constructing(N, RO);
                            if (OF != 16)
                            {
                                StreamReader SR = new StreamReader("orthog" + M.ToString() + ".txt");
                                string all = SR.ReadToEnd();
                                string[] all_with_space = all.Split('\n', ',');
                                Orthogonal_Matrix = new double[M, M];
                                for (int ii = 0; ii < M; ii++)
                                    for (int jj = 0; jj < M; jj++)
                                    {
                                        Orthogonal_Matrix[ii, jj] = Convert.ToDouble(all_with_space[ii + jj]);
                                    }
                            }
                            else
                            {
                                StreamReader SR = new StreamReader("orthog" + N.ToString() + ".txt");
                                string all = SR.ReadToEnd();
                                string[] all_with_space = all.Split('\n', ',');
                                Orthogonal_Matrix = new double[N, N];
                                for (int ii = 0; ii < N; ii++)
                                    for (int jj = 0; jj < N; jj++)
                                    {
                                        Orthogonal_Matrix[ii, jj] = Convert.ToDouble(all_with_space[ii + jj]);
                                    }
                            }
                            //int[] iter = { 10, 20, 30, 40, 50, 100, 200, 300, 500,1000,5000, 10000, 20000 };
                            int[] iter = { 10, 20, 30, 40, 50, 100, 200, 300, 500, 1000, 5000 };
                            //int[] iter = { 1, 2, 3, 4, 5, 6, 7, 8, 9};
                            //10-10 0.255 0.25 0.2460
                            double[] thresolds = { 9000, -500, -5, -43, 30, -2500, 7.5, 65, -4500, -50, -9, -1000, -300, 0, -0.00005, -0.0000001, -0.005, -1600, 44, -1100, -30000, -100000 };
                            string ss = "";
                            string ss2 = "";
                            double[] costs = new double[Ntest];
                            double thresold = thresolds[OF - 1];
                            Random rr = new Random();
                            double Fmax = double.NegativeInfinity; double Fmin = double.PositiveInfinity;
                            for (int IT = 0; IT < iter.Length; IT++)
                            {
                                int[] Fail_Accept = new int[Ntest];
                                int[] iteration_out = new int[Ntest];
                                Iteration = iter[IT];
                                Random Rand = new Random();
                                int FA = -1;
                                int Iter = 0;
                                #region MOA
                                if (method == "MOA")
                                {
                                    double short_range = 0.01;
                                    int acceleration = Convert.ToInt32(Ac.Text);
                                    int Distance = Convert.ToInt32(D.Text);
                                    double[] Magnet = new double[NPop];
                                    double[] Mass = new double[NPop];
                                    double[] Intinsities = { 0.0010, 0.0100, 0.0100, 0.0100, 0.0100, 0.0100, 10.0000, 0.0100, 0.0010, 0.0010, 0.0100, 0.0010, 0.0100, 0, 0.0100, 0.0010, 0.0010, 0.0100, 0.0010, 0.0010, 0.0100, 0.0010 };
                                    double[] Alphas = { 0.0010, 1.0000, 0.0100, 1.0000, 1.0000, 1.0000, 0.0010, 1.0000, 1.0000, 1.0000, 0.0100, 1.0000, 1.0000, 0, 1.0000, 1.0000, 1.0000, 1.0000, 0.0010, 1.0000, 1.0000, 1.0000 };
                                    //for (int alpha_index = 0; alpha_index < 5; alpha_index ++)
                                    //    for (int intensity_index = 0; intensity_index < 5; intensity_index ++)
                                    //    {
                                    double intensity = Intinsities[OF - 1];
                                    double alpha = Alphas[OF - 1];

                                    ss = "Fault_N_" + N.ToString() + "_Pop_" + NPop.ToString() + method + Distance.ToString() + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_OF" + OF.ToString() + "_Th_" + thresold.ToString() + ".txt";
                                    ss2 = "fit_N_" + N.ToString() + "_Pop_" + NPop.ToString() + method + Distance.ToString() + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_OF" + OF.ToString() + "_Th_" + thresold.ToString() + ".txt";
                                    string ss3 = "iter_N_" + N.ToString() + "_Pop_" + NPop.ToString() + method + Distance.ToString() + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_OF" + OF.ToString() + "_Th_" + thresold.ToString() + ".txt";
                                    if (!File.Exists(ss2))
                                    {
                                        for (int tt = 0; tt < Ntest; tt++)
                                        {
                                            FA = -1;
                                            NF[] Nodes = new NF[NPop];
                                            NF[] Forces = new NF[NPop];
                                            NF[] a = new NF[NPop];
                                            NF[] v = new NF[NPop];
                                            NF Best_indi = new NF();
                                            double[] distances = new double[NPop];
                                            Best_indi.F = -1E280;
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                Forces[ii].x = new double[N];
                                                a[ii].x = new double[N];
                                                v[ii].x = new double[N];
                                            }
                                            for (int ii = 0; ii < NPop; ii++)
                                            {
                                                Nodes[ii].x = NFC.node_constructing(N, Rand);
                                            }
                                            for (int it = 0; it < Iteration; it++)
                                            {
                                                Fmax = -1E280; Fmin = 100000000;
                                                int min_index = 0; int max_index = 0;
                                                for (int ii = 0; ii < NPop; ii++)
                                                {
                                                    if (OF < 16 || OF == 22)
                                                        Nodes[ii].F = NFC.compute_fitness(Nodes[ii].x, OF, N);
                                                    else
                                                        Nodes[ii].F = NFC.compute_fitness(Nodes[ii].x, OF, N, Perm, Orthogonal_Matrix, M, o);
                                                    if (Fmax < Nodes[ii].F)
                                                    {
                                                        Fmax = Nodes[ii].F;
                                                        max_index = ii;
                                                    }
                                                    if (Fmin > Nodes[ii].F)
                                                    {
                                                        Fmin = Nodes[ii].F;
                                                        min_index = ii;
                                                    }
                                                }

                                                if (Best_indi.F < Fmax)
                                                {
                                                    Best_indi.F = Fmax;
                                                    Best_indi.x = (double[])Nodes[max_index].x.Clone();
                                                }
                                                double range = Fmax - Fmin;
                                                for (int ii = 0; ii < NPop; ii++)
                                                {
                                                    Magnet[ii] = (Nodes[ii].F - Fmin) / (range);
                                                    Mass[ii] = Magnet[ii] * intensity + alpha;
                                                }
                                                for (int ii = 0; ii < NPop; ii++)
                                                {
                                                    Forces[ii].x = new double[N];
                                                }
                                                for (int ii = 0; ii < NPop; ii++)
                                                {
                                                    int[] Neigbour = NFC.Structures(ii, NPop, structure);
                                                    for (int jj = 0; jj < Neigbour.Length; jj++)
                                                    {
                                                        if (ii != Neigbour[jj])
                                                        {
                                                            if (D.Text == "1")
                                                                distances[ii] = distance1(Nodes[ii].x, Nodes[Neigbour[jj]].x, N);
                                                            //else
                                                            //    distances[ii] = distance4(LHC[ii].LH, LHC[Neigbour[jj]].LH, N, P);
                                                            if (distances[ii] > short_range)
                                                            {
                                                                for (int n = 0; n < N; n++)
                                                                    Forces[ii].x[n] = Forces[ii].x[n] + (Nodes[Neigbour[jj]].x[n] - Nodes[ii].x[n]) * Magnet[Neigbour[jj]] / distances[ii];

                                                            }

                                                            //else
                                                            //{
                                                            //    for (int n = 0; n < N; n++)
                                                            //            Forces[ii].x[n] = rr.NextDouble();

                                                            //    Nodes[ii].x = NFC.node_constructing(N, P, rr);
                                                            //}
                                                        }
                                                    }
                                                    //for (int n = 0; n < N; n++)
                                                    //    for (int p = 0; p < P; p++)
                                                    //    {
                                                    //        Forces[ii].x[n, p] = Forces[ii].x[n, p] / (NPop - 1);
                                                    //        //Forces[ii].x[n, p] = Forces[ii].x[n, p] * rr.NextDouble();
                                                    //    }

                                                }
                                                if (acceleration == 0)
                                                {
                                                    for (int ii = 0; ii < NPop; ii++)
                                                    {
                                                        for (int n = 0; n < N; n++)
                                                        {
                                                            v[ii].x[n] = (Forces[ii].x[n] / Mass[ii]);
                                                            Nodes[ii].x[n] = Nodes[ii].x[n] + v[ii].x[n];
                                                        }

                                                        Nodes[ii].x = NFC.Node_treatment(Nodes[ii].x, N, Rand);
                                                    }
                                                }
                                                else if (acceleration == 1)
                                                {
                                                    for (int ii = 0; ii < NPop; ii++)
                                                    {
                                                        for (int n = 0; n < N; n++)
                                                        {
                                                            a[ii].x[n] = (Forces[ii].x[n] / Mass[ii]);
                                                            v[ii].x[n] = v[ii].x[n] * Rand.NextDouble() + a[ii].x[n];
                                                            Nodes[ii].x[n] = Nodes[ii].x[n] + v[ii].x[n];
                                                        }

                                                        Nodes[ii].x = NFC.Node_treatment(Nodes[ii].x, N, Rand);
                                                    }
                                                }
                                                if (Best_indi.F > thresold)
                                                {
                                                    FA = 1;
                                                    break;
                                                }
                                                Iter = it;
                                            }

                                            costs[tt] = Best_indi.F;
                                            Fail_Accept[tt] = FA;
                                            iteration_out[tt] = (Iter + 1) * NPop;

                                        }

                                        sw = new StreamWriter(ss);
                                        for (int tt = 0; tt < Ntest; tt++)
                                        {
                                            sw.Write(Fail_Accept[tt].ToString() + " ");
                                        }
                                        sw.Close();
                                        sw = new StreamWriter(ss2);
                                        for (int tt = 0; tt < Ntest; tt++)
                                        {
                                            sw.Write(costs[tt].ToString() + " ");
                                        }
                                        sw.Close();
                                        sw = new StreamWriter(ss3);
                                        for (int tt = 0; tt < Ntest; tt++)
                                        {
                                            sw.Write(iteration_out[tt].ToString() + " ");
                                        }
                                        sw.Close();
                                    }
                                    //}

                                }
                                #endregion
                            }
                        }
                    }
                #endregion
            }
            Application.Exit();
        }
        public double distance1(double[] x1, double[] x2, int N)
        {
            double result = 0;
            for (int ii = 0; ii < N; ii++)
                result += Math.Abs(x1[ii] - x2[ii]);
            return result;
        }
        public double distance2(double[,] x1, double[,] x2, int N, int P)
        {
            double result = 0;
            for (int ii = 0; ii < N; ii++)
                for (int jj = 0; jj < P; jj++)
                {
                    result += Math.Pow(Math.Abs(x1[ii, jj] - x2[ii, jj]), 2);
                }
            return Math.Sqrt(result);
        }
        public double distance3(double[,] x1, double[,] x2, int N, int P)
        {
            double result = 0;
            for (int ii = 0; ii < N; ii++)
                for (int jj = 0; jj < P; jj++)
                {
                    result = Math.Max(x1[ii, jj] - x2[ii, jj], result);
                }
            return result;
        }

        public int distance4(int[,] LH1, int[,] LH2, int N, int P)
        {
            int dist = 0;
            for (int ii = 0; ii < N; ii++)
                for (int jj = 0; jj < P; jj++)
                    if (LH1[ii, jj] != LH2[ii, jj])
                        dist++;
            return dist / (P * N);
        }
        public int[,] string_to_array(string ss, int N, int P)
        {
            string[] all = ss.Split(' ');
            int[,] LH = new int[N, P];
            int kk = 0;
            for (int ii = 0; ii < N; ii++)
                for (int jj = 0; jj < P; jj++)
                    LH[ii, jj] = Convert.ToInt32(all[kk++]);
            return LH;
        }
        private void comboBox1_SelectedIndexChanged(object sender, EventArgs e)
        {
            if (comboBox1.Text == "MOA" || comboBox1.Text == "FMOA" || comboBox1.Text == "MA-SW" || comboBox1.Text == "3TMA" || comboBox1.Text == "3SOME" || comboBox1.Text == "BBO" || comboBox1.Text == "JADE" || comboBox1.Text == "PMS" || comboBox1.Text == "MASSW" || comboBox1.Text == "CCPSO2" || comboBox1.Text == "(1+1)LMA" || comboBox1.Text == "2LMA" || comboBox1.Text == "12LMA" || comboBox1.Text == "3LMA_test" || comboBox1.Text == "1LMA" || comboBox1.Text == "MASSW_S" || comboBox1.Text == "FRFMOA" || comboBox1.Text == "OMOA" || comboBox1.Text == "FMOAT"||comboBox1.Text=="CMAES")
            {
                MOA.Enabled = true;
                button1.Enabled = true;
                SA.Enabled = false;
                GA.Enabled = false;
                PSO.Enabled = false;
            }
            else if (comboBox1.Text == "Genetic")
            {
                GA.Enabled = true;
                button1.Enabled = true;
                SA.Enabled = false;
                PSO.Enabled = false;
                MOA.Enabled = false;
            }
            else if (comboBox1.Text == "PSO" || comboBox1.Text == "MAPSO")
            {

                PSO.Enabled = true;
                button1.Enabled = true;
                SA.Enabled = false;
                GA.Enabled = false;
                MOA.Enabled = false;
            }
            else if (comboBox1.Text == "ODE" || comboBox1.Text == "DE")
            {
                ODE.Enabled = true;
                button1.Enabled = true;
            }
            else if (comboBox1.Text == "ES" || comboBox1.Text == "EP" || comboBox1.Text == "FES" || comboBox1.Text == "FEP" || comboBox1.Text == "MA-SW" || comboBox1.Text == "ACO" || comboBox1.Text == "SRE")
            {
                button1.Enabled = true;
            }
        }

        private void button2_Click(object sender, EventArgs e)
        {
            string method = "";
            int Ntest = Convert.ToInt32(NT.Text);
            int NPop = Convert.ToInt32(NP.Text);
            int Iteration = Convert.ToInt32(NI.Text);
            if (comboBox1.Text == "MA-SW")
                method = "MA-SW";
            else if (comboBox1.Text == "3TMA")
                method = "3TMA";
            else if (comboBox1.Text == "3SOME")
                method = "3SOME";
            else if (comboBox1.Text == "BBO")
                method = "BBO";
            else if (comboBox1.Text == "JADE")
                method = "JADE";
            else if (comboBox1.Text == "PMS")
                method = "PMS";
            else if (comboBox1.Text == "MASSW")
                method = "MASSW";
            else if (comboBox1.Text == "CCPSO2")
                method = "CCPSO2";
            else if (comboBox1.Text == "FES")
                method = "FES";
            else if (comboBox1.Text == "PSO")
                method = "PSO";
            else if (comboBox1.Text == "Genetic")
                method = "Genetic";
            else if (comboBox1.Text == "MASSW_S")
                method = "MASSW_S";
            bool nonUniform = false;
            Numerical_functions NFC = new Numerical_functions();
            Random r_constant = new Random(10);
            Random Rand = new Random();
            int No_anchor = 20;
            int No_nonanchor = 180;
            int Area_size = 1;
            double range = Convert.ToDouble(Range.Text);
            double Noise_Factor = 0.01;
            indi ALL = new indi();
            double[] Anchorx = new double[No_anchor];
            double[] Anchory = new double[No_anchor];
            double[] Realx = new double[No_nonanchor];
            double[] Realy = new double[No_nonanchor];
            if (nonUniform)
            {
                ALL = NFC.construct_Non_Anchor_node(No_anchor, r_constant, Area_size);
                Anchorx = (double[])ALL.x.Clone();
                Anchory = (double[])ALL.y.Clone();
                indi NAn = NFC.construct_Non_Anchor_node(No_nonanchor, r_constant, Area_size);
                Realx = (double[])NAn.x.Clone();
                Realy = (double[])NAn.y.Clone();
            }
            else
            {
                NFC.Uniformly_distribution(ref Realx, ref Realy, r_constant, No_nonanchor, No_anchor, Area_size);
                NFC.Uniformly_distribution_anchor_nodes(ref Anchorx, ref Anchory, r_constant, No_nonanchor, No_anchor, Area_size);
            }
            double[] gaussian = new double[No_nonanchor];
            //StreamWriter SW = new StreamWriter("f:\\nodesx.txt");
            //StreamWriter SW2 = new StreamWriter("f:\\nodesy.txt");
            //for (int ii = 0; ii < No_anchor; ii++)
            //{
            //    SW.Write(Anchorx[ii]+" ");
            //    SW2.Write(Anchory[ii] + " ");
            //}
            //for (int jj = 0; jj < No_nonanchor; jj++)
            //{
            //    SW.Write(Realx[jj] + " ");
            //    SW2.Write(Realy[jj] + " ");
            //}
            //SW.Close();
            //SW2.Close();
            //int[] NonAnchAchor = new int[No_nonanchor];
            //for (int ii = 0; ii < No_nonanchor; ii++)
            //{
            //   for(int jj=0;jj<No_anchor;jj++)
            //    {
            //        double range_temp = NFC.compute_distance(Anchorx[jj], Anchory[jj], Realx[ii], Realy[ii]);
            //        if (range_temp < range)
            //        {
            //            NonAnchAchor[ii]++;
            //        }
            //    }
            //}
            //StreamWriter SW = new StreamWriter("f:\\nodes.txt");
            //for (int jj = 0; jj < No_nonanchor; jj++)
            //{
            //    SW.Write(NonAnchAchor[jj] + " ");
            //}
            //SW.Close();

            double[] costs = new double[Ntest];
            for (int ii = 0; ii < gaussian.Length; ii++)
                gaussian[ii] = NFC.Normal_Distribution(r_constant, 0, 1);
            #region 3LMA
            if (method == "3TMA")
            {
                double O = NPop * Convert.ToInt32(NI.Text);
                //double[] N_LS_GL = { 0.3,0.1,0.3,0.9,0.5,0.5,0.7,0.3,0.9,0.1,0.3,0.9,0.7,0,0.9,0.3,0.7,0.5,0.9,0.3,0.9,0.3};
                //double[] Mu_rate_s = { 0.1,0.8,0.5,0.8,0.8,1,1,1,0.8,1,0.8,0.8,0.8,0,1,1,0.8,0.8,1,0.5,1,0.01};
                double[] N_LS_GL = { 0.1, 0.3, 0.5, 0.7, 0.9 };
                //double[] Mu_rate_s = { 0.1, 0.5,0.8,1 };
                //double[] Jrs = { 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };
                //string[] Structs = { "Ring", "Star", "Btree", "Bipartie", "Cluster", "Grid", "Ladder", "Cross-Ladder", "Cellular" };
                string[] Structs = { "Star" };
                //for (int Portion_index = 3; Portion_index < 5; Portion_index++)
                ////    for (int M_index = 0; M_index < 15; M_index++)
                //    {
                indi BEST = new indi();
                bool OBL = true;
                //for (int OB = 0; OB < 10; OB++)
                //{
                for (int SS = 0; SS < 1; SS++)
                {
                    string Struct = Structs[SS];
                    //double[] N_LS_GL = { 0.3, 0.5, 0.9, 0.5, 0.5, 0.3, 0.3, 0.5, 0.5, 0.3, 0.7, 0.5, 0.5, 0, 0.5, 0.7, 0.9, 0.9, 0.9, 0.7, 0.5, 0.7 };
                    //double[] Mu_rate_s = { 0.1, 0.00001, 5, 0.1, 0.1, 0.5, 0.01, 0.1, 5, 0.5, 0.001, 5, 0.1, 0, 0.1, 5, 0.01, 0.5, 0.01, 0.1, 0.5, 0.5 };
                    double LS_LG = 0.5;
                    double Mu_rate = 0.5;
                    Iteration = Convert.ToInt32((O * (1 - LS_LG)) / NPop) + 1;
                    int I_str = Convert.ToInt32(O * LS_LG) / Iteration;
                    int best_index = 0;
                    bool local_search = false;
                    double Jr = 0.3;
                    string ss2 = "fit_Uniform" + "No_A" + No_nonanchor.ToString() + "range" + range.ToString() + "_Mr_" + Mu_rate.ToString() + "LSLG" + LS_LG.ToString() + "_J" + Jr.ToString() + "_Pop_" + NPop.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "Localization" + ".txt";
                    string ss3 = "loc_Uniform" + "No_A" + No_nonanchor.ToString() + "range" + range.ToString() + "_Mr_" + Mu_rate.ToString() + "LSLG" + LS_LG.ToString() + "_J" + Jr.ToString() + "_Pop_" + NPop.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "Localization" + ".txt";
                    if (!File.Exists(ss2))
                    {
                        for (int tt = 0; tt < Ntest; tt++)
                        {

                            //counter_JR = 1;
                            indi[] Pop = new indi[NPop];
                            indi[] OPop = new indi[NPop];
                            indi Best_nodes = new indi();
                            indi[] Breed = new indi[NPop];
                            Best_nodes.cost = 1E250;
                            for (int ii = 0; ii < NPop; ii++)
                            {
                                Pop[ii] = NFC.Initialization(No_nonanchor, Rand);
                                if (OBL)
                                {
                                    OPop[ii].x = NFC.Opposition_Learning(Pop[ii].x);
                                    OPop[ii].y = NFC.Opposition_Learning(Pop[ii].y);
                                }
                                Pop[ii].cost = NFC.compute_cost(Pop[ii].x, Pop[ii].y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);
                                if (OBL)
                                {
                                    OPop[ii].cost = NFC.compute_cost(OPop[ii].x, OPop[ii].y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);

                                    if (OPop[ii].cost < Pop[ii].cost)
                                    {
                                        Pop[ii].cost = OPop[ii].cost;
                                        Pop[ii].x = (double[])OPop[ii].x.Clone();
                                        Pop[ii].y = (double[])OPop[ii].y.Clone();
                                    }
                                }
                                Pop[ii].local_search = true;

                            }
                            for (int it = 0; it < Iteration; it++)
                            {
                                //perform Steady State GA
                                #region generate Breeds
                                double min = 1E250; int min_index = 0;
                                local_search = false;
                                for (int ii = 0; ii < NPop; ii++)
                                {
                                    int parent1 = Rand.Next(NPop);
                                    int[] S = NFC.Structures(ii, NPop, Struct);
                                    int index2 = Rand.Next(S.Length);
                                    //int parent2 = S[index2];
                                    int parent2 = NFC.negative_assortative_mating_strategy(Pop[parent1].x, Pop[parent1].y, Pop, Rand);
                                    Breed[ii].x = NFC.BLX_Crossover(Pop[parent1].x, Pop[parent2].x);
                                    Breed[ii].y = NFC.BLX_Crossover(Pop[parent1].y, Pop[parent2].y);
                                    if (Rand.NextDouble() < Mu_rate)
                                        NFC.Mutation_BGA2(ref Breed[ii].x, ref Breed[ii].y, Pop, No_nonanchor, Rand);
                                    Breed[ii].cost = NFC.compute_cost(Pop[ii].x, Pop[ii].y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);

                                    NFC.Node_treatment_localization(ref Breed[ii].x, ref  Breed[ii].y, No_nonanchor, Rand);
                                    ////////////////////////////////////////////
                                    if (Pop[ii].cost > Breed[ii].cost)
                                    {
                                        Pop[ii].cost = Breed[ii].cost;
                                        Pop[ii].x = (double[])Breed[ii].x.Clone();
                                        Pop[ii].y = (double[])Breed[ii].y.Clone();
                                    }
                                }
                                for (int ii = 0; ii < NPop; ii++)
                                {
                                    if (OBL)
                                    {
                                        if (Jr > Rand.NextDouble())
                                        {

                                            OPop[ii] = NFC.generating_jumping(Pop, Pop[ii].x, NPop, No_nonanchor, ii);
                                            OPop[ii].cost = NFC.compute_cost(Pop[ii].x, Pop[ii].y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);

                                            if (OPop[ii].cost < Pop[ii].cost)
                                            {
                                                Pop[ii].cost = OPop[ii].cost;
                                                Pop[ii].x = (double[])OPop[ii].x.Clone();
                                                Pop[ii].y = (double[])OPop[ii].y.Clone();
                                            }
                                        }
                                    }
                                    if (Pop[ii].cost < min && Pop[ii].local_search == true)
                                    {
                                        min = Pop[ii].cost;
                                        min_index = ii;
                                        local_search = true;
                                    }
                                }
                                #endregion
                                // Replacement Strategy

                                //
                                //Pick the best individual of Population
                                if (local_search == false)
                                {

                                    for (int ii = 0; ii < NPop; ii++)
                                    {
                                        Pop[ii] = NFC.Initialization(No_nonanchor, Rand);
                                        Pop[ii].cost = NFC.compute_cost(Pop[ii].x, Pop[ii].y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);
                                        Pop[ii].local_search = true;
                                        Pop[ii].visited = 0;
                                        if (Pop[ii].cost < min && Pop[ii].local_search == true)
                                        {
                                            min = Pop[ii].cost;
                                            min_index = ii;
                                        }
                                    }
                                    int randx = Rand.Next(0, NPop);
                                    Pop[randx].x = (double[])Best_nodes.x.Clone();
                                    Pop[randx].y = (double[])Best_nodes.y.Clone();
                                    Pop[randx].cost = Best_nodes.cost;
                                }
                                indi Best_Indi_BLS = new indi();
                                Best_Indi_BLS.cost = Pop[min_index].cost;
                                Best_Indi_BLS.x = (double[])Pop[min_index].x.Clone();
                                Best_Indi_BLS.y = (double[])Pop[min_index].y.Clone();
                                Best_Indi_BLS.visited = Pop[min_index].visited;
                                Best_Indi_BLS.rho = Pop[min_index].rho;
                                if (Pop[min_index].visited == 1)
                                {
                                    Best_Indi_BLS.biasx = (double[])Pop[min_index].biasx.Clone();
                                    Best_Indi_BLS.biasy = (double[])Pop[min_index].biasy.Clone();
                                }
                                best_index = min_index;
                                // Initialize Bias and Rho
                                if (Best_Indi_BLS.visited == 0)
                                {
                                    Best_Indi_BLS.biasx = new double[No_nonanchor];
                                    Best_Indi_BLS.biasy = new double[No_nonanchor];
                                    Best_Indi_BLS.rho = 0.01;
                                }
                                // Local Search Operator
                                indi Best_Indi_ALS = new indi();
                                Best_Indi_ALS.x = (double[])Best_Indi_BLS.x.Clone();
                                Best_Indi_ALS.y = (double[])Best_Indi_BLS.y.Clone();
                                Best_Indi_ALS.cost = Best_Indi_BLS.cost;
                                Best_Indi_ALS.rho = Best_Indi_BLS.rho;
                                Best_Indi_ALS.visited = Best_Indi_BLS.visited;
                                Best_Indi_ALS.biasx = (double[])Best_Indi_BLS.biasx.Clone();
                                Best_Indi_ALS.biasy = (double[])Best_Indi_BLS.biasy.Clone();
                                NFC.SolisWets(ref Best_Indi_ALS.x, ref Best_Indi_ALS.y, ref Best_Indi_ALS.cost, No_nonanchor, ref Best_Indi_ALS.biasx, ref Best_Indi_ALS.biasy, ref Best_Indi_ALS.rho, I_str, Rand, range, Realx, Realy, gaussian, Noise_Factor, Anchorx, Anchory);
                                if (Best_Indi_BLS.cost - Best_Indi_ALS.cost == 0)
                                {
                                    Pop[best_index].local_search = false;
                                    //  Pop[best_index].visited = 1;
                                }
                                else
                                {
                                    Pop[best_index].x = (double[])Best_Indi_ALS.x.Clone();
                                    Pop[best_index].y = (double[])Best_Indi_ALS.y.Clone();
                                    Pop[best_index].visited = 1;
                                    Pop[best_index].cost = Best_Indi_ALS.cost;
                                    Pop[best_index].rho = Best_Indi_ALS.rho;
                                    Pop[best_index].biasx = (double[])Best_Indi_ALS.biasx.Clone();
                                    Pop[best_index].biasy = (double[])Best_Indi_ALS.biasy.Clone();
                                    Pop[best_index].local_search = true;
                                }

                                if (Best_nodes.cost > Pop[best_index].cost)
                                {
                                    Best_nodes.cost = Pop[best_index].cost;
                                    Best_nodes.x = (double[])Pop[best_index].x.Clone();
                                    Best_nodes.y = (double[])Pop[best_index].y.Clone();
                                }

                                /////////////////////
                                //counter_JR=counter_JR * 0.99;
                                //Jr = 1 - (counter_JR);
                            }
                            costs[tt] = NFC.Calc_error(Best_nodes.x, Best_nodes.y, Realx, Realy, range);
                            BEST.x = (double[])Best_nodes.x.Clone();
                            BEST.y = (double[])Best_nodes.y.Clone();
                        }

                        StreamWriter sw2 = new StreamWriter(ss3);

                        for (int ii = 0; ii < No_nonanchor; ii++)
                        {
                            sw2.Write(BEST.x[ii] + " ");
                        }
                        sw2.WriteLine();
                        for (int ii = 0; ii < No_nonanchor; ii++)
                        {
                            sw2.Write(BEST.y[ii] + " ");
                        }
                        sw2.Close();
                        sw = new StreamWriter(ss2);
                        for (int tt = 0; tt < Ntest; tt++)
                        {
                            sw.Write(costs[tt].ToString() + " ");
                        }
                        sw.Close();
                    }
                }
                //}
            }

            #endregion
            #region 3SOME
            if (method == "3SOME")
            {
                indi BEST = new indi();
                int k = 4;
                double[] inheritance_factors = { 0.3, 0.6, 0.3, 0.6, 0.6, 0.9, 0.3, 0.6, 0.6, 0.3, 0.6, 0.9, 0.9, 0, 0.6, 0.9, 0.6, 0.1, 0.3, 0.9, 0.01, 0.6 };
                //double[] deltas = { 0.001, 0.01, 0.4, 0.001, 0.001, 0.01, 0.4, 0.01, 0.001, 0.1, 0.01, 0.01, 0.001, 0, 0.001, 0.001, 0.01, 0.01, 0.001, 0.1, 0.01, 0.01 };
                //double[] rhos = { 0.01, 0.1, 0.2, 0.4, 0.6 };
                bool updated = false;
                //for (int inh = 0; inh < 5; inh++)
                double[] deltas = { 0.001, 0.01 };
                for (int del = 0; del < 1; del++)
                {

                    double rho = 0.4;
                    double inheritance_factor = 0.6;
                    double delta = deltas[del];
                    double root = No_nonanchor * inheritance_factor;
                    double Cr = 1.0 / (Math.Pow(2, 1.0 / root));
                    root = No_nonanchor * (1 - inheritance_factor);
                    double Crprim = 1.0 / (Math.Pow(2, 1.0 / root));
                    int local_budget = 150;
                    string ss2 = "fit_N_Uniform" + No_nonanchor.ToString() + "pop_" + NPop.ToString() + "Local" + method + "_inh_" + inheritance_factor.ToString() + "_del_" + delta.ToString() + "_rho_" + rho.ToString() + "_test_" + Ntest.ToString() + "range" + range.ToString() + ".txt";
                    string ss3 = "loc_N_Uniform" + No_nonanchor.ToString() + "pop_" + NPop.ToString() + "Local" + method + "_inh_" + inheritance_factor.ToString() + "_del_" + delta.ToString() + "_rho_" + rho.ToString() + "_test_" + Ntest.ToString() + "range" + range.ToString() + ".txt";
                    #region file_existency
                    if (!File.Exists(ss2))
                    {
                        for (int tt = 0; tt < Ntest; tt++)
                        {
                            rho = 0.4;
                            indi Elite = new indi();
                            indi Indi = new indi();
                            ////////////////////////////////////
                            //////////////////////////////
                            Elite = NFC.Initialization(No_nonanchor, Rand);
                            Elite.cost = NFC.compute_cost(Elite.x, Elite.y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);
                            for (int it = 0; it < Iteration * NPop; )
                            {
                                ///////////////////Long distance/////////////////////////////////
                                updated = false;
                                while (!updated && it < Iteration * NPop)
                                {
                                    #region Long_distance
                                    int checker = 0;
                                    Indi = NFC.Initialization(No_nonanchor, Rand);
                                    int i = Rand.Next(No_nonanchor);
                                    Indi.x[i] = Elite.x[i]; Indi.y[i] = Elite.y[i];
                                    while (Rand.NextDouble() <= Cr)
                                    {
                                        Indi.x[i] = Elite.x[i];
                                        Indi.y[i] = Elite.y[i];
                                        i++;
                                        if (i == No_nonanchor)
                                            i = 0;
                                        if (checker >= No_nonanchor)
                                            break;
                                        checker++;
                                    }
                                    Indi.cost = NFC.compute_cost(Indi.x, Indi.y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);
                                    if (Indi.cost <= Elite.cost)
                                    {
                                        Elite.x = (double[])Indi.x.Clone();
                                        Elite.y = (double[])Indi.y.Clone();
                                        Elite.cost = Indi.cost;
                                        updated = true;
                                    }
                                    //////////////////////Function Evaluation increasing///////////////////
                                    it++;
                                    #endregion
                                }
                                ///////////////Middle Distance////////////////////////////
                                while (updated && it < Iteration * NPop)
                                {
                                    updated = false;
                                    #region Middle_distance
                                    for (int I = 0; I < k * No_nonanchor; I++)
                                    {
                                        Indi.x = NFC.Initialization_Hypercube(Elite.x, No_nonanchor, Rand, delta);
                                        Indi.y = NFC.Initialization_Hypercube(Elite.y, No_nonanchor, Rand, delta);
                                        int i = Rand.Next(No_nonanchor);
                                        Indi.x[i] = Elite.x[i];
                                        Indi.y[i] = Elite.y[i];
                                        int checker = 0;
                                        while (Rand.NextDouble() <= Crprim)
                                        {
                                            Indi.x[i] = Elite.x[i];
                                            Indi.y[i] = Elite.y[i];
                                            i++;
                                            if (i == No_nonanchor)
                                                i = 0;
                                            if (checker >= No_nonanchor)
                                                break;
                                            checker++;
                                        }
                                        NFC.Node_treatment_localization(ref Indi.x, ref Indi.y, No_nonanchor, Rand);
                                        Indi.cost = NFC.compute_cost(Indi.x, Indi.y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);
                                        if (Indi.cost < Elite.cost)
                                        {
                                            updated = true;
                                        }
                                        if (Indi.cost <= Elite.cost)
                                        {
                                            Elite.x = (double[])Indi.x.Clone();
                                            Elite.cost = Indi.cost;
                                        }
                                        it++;
                                    }
                                    #endregion
                                }
                                /////////////////////////Short Distance////////////////////////////////////////
                                #region Short_region
                                indi indi2 = new indi();
                                for (int ii = 0; ii < local_budget && it < NPop * Iteration; ii++)
                                {
                                    indi2.x = (double[])Elite.x.Clone();
                                    indi2.y = (double[])Elite.y.Clone();
                                    for (int jj = 0; jj < No_nonanchor; jj++)
                                    {
                                        if (Elite.x[jj] - rho > -1)
                                        {
                                            indi2.x[jj] = Elite.x[jj] - rho;
                                            indi2.y[jj] = Elite.y[jj] - rho;
                                        }
                                        else
                                        {
                                            indi2.x[jj] = Elite.x[jj];
                                            indi2.y[jj] = Elite.y[jj];
                                        }
                                        indi2.cost = NFC.compute_cost(indi2.x, indi2.y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);
                                        it++;
                                        if (indi2.cost <= Indi.cost)
                                        {
                                            Indi.x = (double[])indi2.x.Clone();
                                            Indi.y = (double[])indi2.y.Clone();
                                            Indi.cost = indi2.cost;
                                        }
                                        else
                                        {
                                            if (Elite.x[jj] + rho / 2 < 1)
                                                indi2.x[jj] = Elite.x[jj] + rho / 2;
                                            else
                                                indi2.x[jj] = Elite.x[jj];
                                            indi2.cost = NFC.compute_cost(indi2.x, indi2.y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);
                                            it++;
                                            if (indi2.cost <= Indi.cost)
                                            {
                                                Indi.x = (double[])indi2.x.Clone();
                                                Indi.y = (double[])indi2.y.Clone();
                                                Indi.cost = indi2.cost;
                                            }
                                        }
                                    }
                                    if (Indi.cost <= Elite.cost)
                                    {
                                        Elite.x = (double[])Indi.x.Clone();
                                        Elite.y = (double[])Indi.y.Clone();
                                        Elite.cost = Indi.cost;
                                        updated = true;
                                    }
                                    else
                                        rho /= 2;

                                }
                                #endregion
                                if (it < NPop * Iteration)
                                {
                                    if (updated)
                                    {
                                        #region Middle_distance
                                        for (int I = 0; I < k * No_nonanchor; I++)
                                        {
                                            Indi.x = NFC.Initialization_Hypercube(Elite.x, No_nonanchor, Rand, delta);
                                            Indi.y = NFC.Initialization_Hypercube(Elite.y, No_nonanchor, Rand, delta);
                                            int i = Rand.Next(No_nonanchor);
                                            Indi.x[i] = Elite.x[i];
                                            Indi.y[i] = Elite.y[i];
                                            int checker = 0;
                                            while (Rand.NextDouble() <= Crprim)
                                            {
                                                Indi.x[i] = Elite.x[i];
                                                Indi.y[i] = Elite.y[i];
                                                i++;
                                                if (i == No_nonanchor)
                                                    i = 0;
                                                if (checker >= No_nonanchor)
                                                    break;
                                                checker++;
                                            }
                                            NFC.Node_treatment_localization(ref Indi.x, ref Indi.y, No_nonanchor, Rand);
                                            Indi.cost = NFC.compute_cost(Indi.x, Indi.y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);
                                            if (Indi.cost < Elite.cost)
                                            {
                                                updated = true;
                                            }
                                            if (Indi.cost <= Elite.cost)
                                            {
                                                Elite.x = (double[])Indi.x.Clone();
                                                Elite.cost = Indi.cost;
                                            }
                                            it++;
                                        }
                                        #endregion
                                    }
                                    else
                                    {
                                        #region Long_distance
                                        int checker = 0;
                                        Indi = NFC.Initialization(No_nonanchor, Rand);
                                        int i = Rand.Next(No_nonanchor);
                                        Indi.x[i] = Elite.x[i]; Indi.y[i] = Elite.y[i];
                                        while (Rand.NextDouble() <= Cr)
                                        {
                                            Indi.x[i] = Elite.x[i];
                                            Indi.y[i] = Elite.y[i];
                                            i++;
                                            if (i == No_nonanchor)
                                                i = 0;
                                            if (checker >= No_nonanchor)
                                                break;
                                            checker++;
                                        }
                                        Indi.cost = NFC.compute_cost(Indi.x, Indi.y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);
                                        if (Indi.cost <= Elite.cost)
                                        {
                                            Elite.x = (double[])Indi.x.Clone();
                                            Elite.y = (double[])Indi.y.Clone();
                                            Elite.cost = Indi.cost;
                                            updated = true;
                                        }
                                        //////////////////////Function Evaluation increasing///////////////////
                                        it++;
                                        #endregion
                                    }
                                }
                            }
                            costs[tt] = NFC.Calc_error(Elite.x, Elite.y, Realx, Realy, range);
                            BEST.x = (double[])Elite.x.Clone();
                            BEST.y = (double[])Elite.y.Clone();
                        }

                        StreamWriter sw2 = new StreamWriter(ss3);

                        for (int ii = 0; ii < No_nonanchor; ii++)
                        {
                            sw2.Write(BEST.x[ii] + " ");
                        }
                        sw2.WriteLine();
                        for (int ii = 0; ii < No_nonanchor; ii++)
                        {
                            sw2.Write(BEST.y[ii] + " ");
                        }
                        sw2.Close();
                        sw = new StreamWriter(ss2);
                        for (int tt = 0; tt < Ntest; tt++)
                        {
                            sw.Write(costs[tt].ToString() + " ");
                        }
                        sw.Close();

                    }
                }
                    #endregion


            }
            #endregion
            #region MASW
            if (method == "MA-SW")
            {
                double O = NPop * Convert.ToInt32(NI.Text);
                //double[] N_LS_GL = { 0.3,0.1,0.3,0.9,0.5,0.5,0.7,0.3,0.9,0.1,0.3,0.9,0.7,0,0.9,0.3,0.7,0.5,0.9,0.3,0.9,0.3};
                //double[] Mu_rate_s = { 0.1,0.8,0.5,0.8,0.8,1,1,1,0.8,1,0.8,0.8,0.8,0,1,1,0.8,0.8,1,0.5,1,0.01};
                double[] N_LS_GL = { 0.1, 0.3, 0.5, 0.7, 0.9 };
                //double[] Mu_rate_s = { 0.1, 0.5,0.8,1 };
                //double[] Jrs = { 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };
                //string[] Structs = { "Ring", "Star", "Btree", "Bipartie", "Cluster", "Grid", "Ladder", "Cross-Ladder", "Cellular" };
                string[] Structs = { "Star" };
                //for (int Portion_index = 1; Portion_index < 5; Portion_index++)
                ////    for (int M_index = 0; M_index < 15; M_index++)
                //{

                indi BEST = new indi();
                //for (int OB = 0; OB < 10; OB++)
                //{
                for (int SS = 0; SS < 1; SS++)
                {
                    string Struct = Structs[SS];
                    //double[] N_LS_GL = { 0.3, 0.5, 0.9, 0.5, 0.5, 0.3, 0.3, 0.5, 0.5, 0.3, 0.7, 0.5, 0.5, 0, 0.5, 0.7, 0.9, 0.9, 0.9, 0.7, 0.5, 0.7 };
                    //double[] Mu_rate_s = { 0.1, 0.00001, 5, 0.1, 0.1, 0.5, 0.01, 0.1, 5, 0.5, 0.001, 5, 0.1, 0, 0.1, 5, 0.01, 0.5, 0.01, 0.1, 0.5, 0.5 };
                    double LS_LG = 0.7;
                    double Mu_rate = 0.5;
                    Iteration = Convert.ToInt32((O * (1 - LS_LG)) / NPop) + 1;
                    int I_str = Convert.ToInt32(O * LS_LG) / Iteration;
                    int best_index = 0;
                    bool local_search = false;
                    double Jr = 0.3;
                    string ss2 = "fit_" + "No_UniformA" + No_nonanchor.ToString() + "range" + range.ToString() + "_Mr_" + Mu_rate.ToString() + "LSLG" + LS_LG.ToString() + "_J" + Jr.ToString() + "_Pop_" + NPop.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "Localization" + ".txt";
                    string ss3 = "loc_" + "No_UniformA" + No_nonanchor.ToString() + "range" + range.ToString() + "_Mr_" + Mu_rate.ToString() + "LSLG" + LS_LG.ToString() + "_J" + Jr.ToString() + "_Pop_" + NPop.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "Localization" + ".txt";
                    if (!File.Exists(ss2))
                    {
                        for (int tt = 0; tt < Ntest; tt++)
                        {

                            //counter_JR = 1;
                            indi[] Pop = new indi[NPop];
                            indi Best_nodes = new indi();
                            indi[] Breed = new indi[NPop];
                            Best_nodes.cost = 1E250;
                            for (int ii = 0; ii < NPop; ii++)
                            {
                                Pop[ii] = NFC.Initialization(No_nonanchor, Rand);
                                Pop[ii].cost = NFC.compute_cost(Pop[ii].x, Pop[ii].y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);
                                Pop[ii].local_search = true;

                            }
                            for (int it = 0; it < Iteration; it++)
                            {
                                //perform Steady State GA
                                #region generate Breeds
                                double min = 1E250; int min_index = 0;
                                local_search = false;
                                int imin = -1; int imax = -1; double min_b = 1E250; double max_p = -1E250;
                                for (int ii = 0; ii < NPop; ii++)
                                {
                                    int parent1 = Rand.Next(NPop);
                                    int[] S = NFC.Structures(ii, NPop, Struct);
                                    int index2 = Rand.Next(S.Length);
                                    //int parent2 = S[index2];
                                    int parent2 = NFC.negative_assortative_mating_strategy(Pop[parent1].x, Pop[parent1].y, Pop, Rand);
                                    Breed[ii].x = NFC.BLX_Crossover(Pop[parent1].x, Pop[parent2].x);
                                    Breed[ii].y = NFC.BLX_Crossover(Pop[parent1].y, Pop[parent2].y);
                                    if (Rand.NextDouble() < Mu_rate)
                                        NFC.Mutation_BGA2(ref Breed[ii].x, ref Breed[ii].y, Pop, No_nonanchor, Rand);
                                    Breed[ii].cost = NFC.compute_cost(Pop[ii].x, Pop[ii].y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);

                                    NFC.Node_treatment_localization(ref Breed[ii].x, ref  Breed[ii].y, No_nonanchor, Rand);
                                    ////////////////////////////////////////////
                                    if (Breed[ii].cost < min_b)
                                    {
                                        imin = ii;
                                        min_b = Breed[ii].cost;
                                    }
                                    if (Pop[ii].cost > max_p)
                                    {
                                        imax = ii;
                                        max_p = Pop[ii].cost;
                                    }
                                    if (Pop[ii].cost < min && Pop[ii].local_search == true)
                                    {
                                        min = Pop[ii].cost;
                                        min_index = ii;
                                        local_search = true;
                                    }
                                }
                                if (min_b < max_p)
                                {
                                    Pop[imax].cost = Breed[imin].cost;
                                    Pop[imax].x = (double[])Breed[imin].x.Clone();
                                    Pop[imax].y = (double[])Breed[imin].y.Clone();
                                }
                                #endregion
                                // Replacement Strategy

                                //
                                //Pick the best individual of Population
                                if (local_search == false)
                                {

                                    for (int ii = 0; ii < NPop; ii++)
                                    {
                                        Pop[ii] = NFC.Initialization(No_nonanchor, Rand);
                                        Pop[ii].cost = NFC.compute_cost(Pop[ii].x, Pop[ii].y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);
                                        Pop[ii].local_search = true;
                                        Pop[ii].visited = 0;
                                        if (Pop[ii].cost < min && Pop[ii].local_search == true)
                                        {
                                            min = Pop[ii].cost;
                                            min_index = ii;
                                        }
                                    }
                                    int randx = Rand.Next(0, NPop);
                                    Pop[randx].x = (double[])Best_nodes.x.Clone();
                                    Pop[randx].y = (double[])Best_nodes.y.Clone();
                                    Pop[randx].cost = Best_nodes.cost;
                                }
                                indi Best_Indi_BLS = new indi();
                                Best_Indi_BLS.cost = Pop[min_index].cost;
                                Best_Indi_BLS.x = (double[])Pop[min_index].x.Clone();
                                Best_Indi_BLS.y = (double[])Pop[min_index].y.Clone();
                                Best_Indi_BLS.visited = Pop[min_index].visited;
                                Best_Indi_BLS.rho = Pop[min_index].rho;
                                if (Pop[min_index].visited == 1)
                                {
                                    Best_Indi_BLS.biasx = (double[])Pop[min_index].biasx.Clone();
                                    Best_Indi_BLS.biasy = (double[])Pop[min_index].biasy.Clone();
                                }
                                best_index = min_index;
                                // Initialize Bias and Rho
                                if (Best_Indi_BLS.visited == 0)
                                {
                                    Best_Indi_BLS.biasx = new double[No_nonanchor];
                                    Best_Indi_BLS.biasy = new double[No_nonanchor];
                                    Best_Indi_BLS.rho = 0.01;
                                }
                                // Local Search Operator
                                indi Best_Indi_ALS = new indi();
                                Best_Indi_ALS.x = (double[])Best_Indi_BLS.x.Clone();
                                Best_Indi_ALS.y = (double[])Best_Indi_BLS.y.Clone();
                                Best_Indi_ALS.cost = Best_Indi_BLS.cost;
                                Best_Indi_ALS.rho = Best_Indi_BLS.rho;
                                Best_Indi_ALS.visited = Best_Indi_BLS.visited;
                                Best_Indi_ALS.biasx = (double[])Best_Indi_BLS.biasx.Clone();
                                Best_Indi_ALS.biasy = (double[])Best_Indi_BLS.biasy.Clone();
                                NFC.SolisWets(ref Best_Indi_ALS.x, ref Best_Indi_ALS.y, ref Best_Indi_ALS.cost, No_nonanchor, ref Best_Indi_ALS.biasx, ref Best_Indi_ALS.biasy, ref Best_Indi_ALS.rho, I_str, Rand, range, Realx, Realy, gaussian, Noise_Factor, Anchorx, Anchory);
                                if (Best_Indi_BLS.cost - Best_Indi_ALS.cost == 0)
                                {
                                    Pop[best_index].local_search = false;
                                    //  Pop[best_index].visited = 1;
                                }
                                else
                                {
                                    Pop[best_index].x = (double[])Best_Indi_ALS.x.Clone();
                                    Pop[best_index].y = (double[])Best_Indi_ALS.y.Clone();
                                    Pop[best_index].visited = 1;
                                    Pop[best_index].cost = Best_Indi_ALS.cost;
                                    Pop[best_index].rho = Best_Indi_ALS.rho;
                                    Pop[best_index].biasx = (double[])Best_Indi_ALS.biasx.Clone();
                                    Pop[best_index].biasy = (double[])Best_Indi_ALS.biasy.Clone();
                                    Pop[best_index].local_search = true;
                                }

                                if (Best_nodes.cost > Pop[best_index].cost)
                                {
                                    Best_nodes.cost = Pop[best_index].cost;
                                    Best_nodes.x = (double[])Pop[best_index].x.Clone();
                                    Best_nodes.y = (double[])Pop[best_index].y.Clone();
                                }

                                /////////////////////
                                //counter_JR=counter_JR * 0.99;
                                //Jr = 1 - (counter_JR);
                            }
                            costs[tt] = NFC.Calc_error(Best_nodes.x, Best_nodes.y, Realx, Realy, range);
                            BEST.x = (double[])Best_nodes.x.Clone();
                            BEST.y = (double[])Best_nodes.y.Clone();
                        }

                        StreamWriter sw2 = new StreamWriter(ss3);

                        for (int ii = 0; ii < No_nonanchor; ii++)
                        {
                            sw2.Write(BEST.x[ii] + " ");
                        }
                        sw2.WriteLine();
                        for (int ii = 0; ii < No_nonanchor; ii++)
                        {
                            sw2.Write(BEST.y[ii] + " ");
                        }
                        sw2.Close();
                        sw = new StreamWriter(ss2);
                        for (int tt = 0; tt < Ntest; tt++)
                        {
                            sw.Write(costs[tt].ToString() + " ");
                        }
                        sw.Close();
                    }
                }
                //}
            }

            #endregion
            #region MASSW
            if (method == "MASSW")
            {
                double O = NPop * Convert.ToInt32(NI.Text);
                //double[] N_LS_GL = { 0.3,0.1,0.3,0.9,0.5,0.5,0.7,0.3,0.9,0.1,0.3,0.9,0.7,0,0.9,0.3,0.7,0.5,0.9,0.3,0.9,0.3};
                //double[] Mu_rate_s = { 0.1,0.8,0.5,0.8,0.8,1,1,1,0.8,1,0.8,0.8,0.8,0,1,1,0.8,0.8,1,0.5,1,0.01};
                double[] N_LS_GL = { 0.1, 0.3, 0.5, 0.7, 0.9 };
                //double[] Mu_rate_s = { 0.1, 0.5,0.8,1 };
                //double[] Jrs = { 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };
                //string[] Structs = { "Ring", "Star", "Btree", "Bipartie", "Cluster", "Grid", "Ladder", "Cross-Ladder", "Cellular" };
                string[] Structs = { "Star" };
                //for (int Portion_index = 1; Portion_index < 3; Portion_index++)
                ////    for (int M_index = 0; M_index < 15; M_index++)
                //{

                indi BEST = new indi();
                //for (int OB = 0; OB < 10; OB++)
                //{
                for (int SS = 0; SS < 1; SS++)
                {
                    string Struct = Structs[SS];
                    //double[] N_LS_GL = { 0.3, 0.5, 0.9, 0.5, 0.5, 0.3, 0.3, 0.5, 0.5, 0.3, 0.7, 0.5, 0.5, 0, 0.5, 0.7, 0.9, 0.9, 0.9, 0.7, 0.5, 0.7 };
                    //double[] Mu_rate_s = { 0.1, 0.00001, 5, 0.1, 0.1, 0.5, 0.01, 0.1, 5, 0.5, 0.001, 5, 0.1, 0, 0.1, 5, 0.01, 0.5, 0.01, 0.1, 0.5, 0.5 };
                    double LS_LG = 0.3;
                    double Mu_rate = 0.01;
                    Iteration = Convert.ToInt32((O * (1 - LS_LG)) / NPop) + 1;
                    int I_str = Convert.ToInt32(O * LS_LG) / Iteration;
                    int best_index = 0;
                    bool local_search = false;
                    string ss2 = "fit_Uniform" + "No_A" + No_nonanchor.ToString() + "range" + range.ToString() + "_Mr_" + Mu_rate.ToString() + "LSLG" + LS_LG.ToString() + "_Pop_" + NPop.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "Localization" + ".txt";
                    string ss3 = "loc_Uniform" + "No_A" + No_nonanchor.ToString() + "range" + range.ToString() + "_Mr_" + Mu_rate.ToString() + "LSLG" + LS_LG.ToString() + "_Pop_" + NPop.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "Localization" + ".txt";
                    if (!File.Exists(ss2))
                    {
                        for (int tt = 0; tt < Ntest; tt++)
                        {

                            //counter_JR = 1;
                            indi[] Pop = new indi[NPop];
                            indi Best_nodes = new indi();
                            indi[] Breed = new indi[NPop];
                            Best_nodes.cost = 1E250;
                            for (int ii = 0; ii < NPop; ii++)
                            {
                                Pop[ii] = NFC.Initialization(No_nonanchor, Rand);
                                Pop[ii].cost = NFC.compute_cost(Pop[ii].x, Pop[ii].y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);
                                Pop[ii].local_search = true;

                            }
                            for (int it = 0; it < Iteration; it++)
                            {
                                //perform Steady State GA
                                #region generate Breeds
                                double min = 1E250; int min_index = 0;
                                local_search = false;
                                int imin = -1; int imax = -1; double min_b = 1E250; double max_p = -1E250;
                                for (int ii = 0; ii < NPop; ii++)
                                {
                                    int parent1 = Rand.Next(NPop);
                                    int[] S = NFC.Structures(ii, NPop, Struct);
                                    int index2 = Rand.Next(S.Length);
                                    //int parent2 = S[index2];
                                    int parent2 = NFC.negative_assortative_mating_strategy(Pop[parent1].x, Pop[parent1].y, Pop, Rand);
                                    Breed[ii].x = NFC.BLX_Crossover(Pop[parent1].x, Pop[parent2].x);
                                    Breed[ii].y = NFC.BLX_Crossover(Pop[parent1].y, Pop[parent2].y);
                                    if (Rand.NextDouble() < Mu_rate)
                                        NFC.Mutation_BGA2(ref Breed[ii].x, ref Breed[ii].y, Pop, No_nonanchor, Rand);
                                    Breed[ii].cost = NFC.compute_cost(Pop[ii].x, Pop[ii].y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);

                                    NFC.Node_treatment_localization(ref Breed[ii].x, ref  Breed[ii].y, No_nonanchor, Rand);
                                    ////////////////////////////////////////////
                                    if (Breed[ii].cost < min_b)
                                    {
                                        imin = ii;
                                        min_b = Breed[ii].cost;
                                    }
                                    if (Pop[ii].cost > max_p)
                                    {
                                        imax = ii;
                                        max_p = Pop[ii].cost;
                                    }
                                    if (Pop[ii].cost < min && Pop[ii].local_search == true)
                                    {
                                        min = Pop[ii].cost;
                                        min_index = ii;
                                        local_search = true;
                                    }
                                }
                                if (min_b < max_p)
                                {
                                    Pop[imax].cost = Breed[imin].cost;
                                    Pop[imax].x = (double[])Breed[imin].x.Clone();
                                    Pop[imax].y = (double[])Breed[imin].y.Clone();
                                }
                                #endregion
                                // Replacement Strategy

                                //
                                //Pick the best individual of Population
                                if (local_search == false)
                                {

                                    for (int ii = 0; ii < NPop; ii++)
                                    {
                                        Pop[ii] = NFC.Initialization(No_nonanchor, Rand);
                                        Pop[ii].cost = NFC.compute_cost(Pop[ii].x, Pop[ii].y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);
                                        Pop[ii].local_search = true;
                                        Pop[ii].visited = 0;
                                        if (Pop[ii].cost < min && Pop[ii].local_search == true)
                                        {
                                            min = Pop[ii].cost;
                                            min_index = ii;
                                        }
                                    }
                                    int randx = Rand.Next(0, NPop);
                                    Pop[randx].x = (double[])Best_nodes.x.Clone();
                                    Pop[randx].y = (double[])Best_nodes.y.Clone();
                                    Pop[randx].cost = Best_nodes.cost;
                                }
                                indi Best_Indi_BLS = new indi();
                                Best_Indi_BLS.cost = Pop[min_index].cost;
                                Best_Indi_BLS.x = (double[])Pop[min_index].x.Clone();
                                Best_Indi_BLS.y = (double[])Pop[min_index].y.Clone();
                                Best_Indi_BLS.visited = Pop[min_index].visited;
                                Best_Indi_BLS.rho = Pop[min_index].rho;
                                if (Pop[min_index].visited == 1)
                                {
                                    Best_Indi_BLS.biasx = (double[])Pop[min_index].biasx.Clone();
                                    Best_Indi_BLS.biasy = (double[])Pop[min_index].biasy.Clone();
                                }
                                best_index = min_index;
                                // Initialize Bias and Rho
                                if (Best_Indi_BLS.visited == 0)
                                {
                                    Best_Indi_BLS.biasx = new double[No_nonanchor];
                                    Best_Indi_BLS.biasy = new double[No_nonanchor];
                                    Best_Indi_BLS.rho = 0.01;
                                }
                                // Local Search Operator
                                indi Best_Indi_ALS = new indi();
                                Best_Indi_ALS.x = (double[])Best_Indi_BLS.x.Clone();
                                Best_Indi_ALS.y = (double[])Best_Indi_BLS.y.Clone();
                                Best_Indi_ALS.cost = Best_Indi_BLS.cost;
                                Best_Indi_ALS.rho = Best_Indi_BLS.rho;
                                Best_Indi_ALS.visited = Best_Indi_BLS.visited;
                                Best_Indi_ALS.biasx = (double[])Best_Indi_BLS.biasx.Clone();
                                Best_Indi_ALS.biasy = (double[])Best_Indi_BLS.biasy.Clone();
                                NFC.Subgrouping_SolisWets_localization(ref Best_Indi_ALS.x, ref Best_Indi_ALS.y, ref Best_Indi_ALS.cost, No_nonanchor, ref Best_Indi_ALS.biasx, ref Best_Indi_ALS.biasy, ref Best_Indi_ALS.rho, I_str, Rand, range, Realx, Realy, gaussian, Noise_Factor, Anchorx, Anchory);
                                if (Best_Indi_BLS.cost - Best_Indi_ALS.cost == 0)
                                {
                                    Pop[best_index].local_search = false;
                                    //  Pop[best_index].visited = 1;
                                }
                                else
                                {
                                    Pop[best_index].x = (double[])Best_Indi_ALS.x.Clone();
                                    Pop[best_index].y = (double[])Best_Indi_ALS.y.Clone();
                                    Pop[best_index].visited = 1;
                                    Pop[best_index].cost = Best_Indi_ALS.cost;
                                    Pop[best_index].rho = Best_Indi_ALS.rho;
                                    Pop[best_index].biasx = (double[])Best_Indi_ALS.biasx.Clone();
                                    Pop[best_index].biasy = (double[])Best_Indi_ALS.biasy.Clone();
                                    Pop[best_index].local_search = true;
                                }

                                if (Best_nodes.cost > Pop[best_index].cost)
                                {
                                    Best_nodes.cost = Pop[best_index].cost;
                                    Best_nodes.x = (double[])Pop[best_index].x.Clone();
                                    Best_nodes.y = (double[])Pop[best_index].y.Clone();
                                }
                            }
                            costs[tt] = NFC.Calc_error(Best_nodes.x, Best_nodes.y, Realx, Realy, range);
                            BEST.x = (double[])Best_nodes.x.Clone();
                            BEST.y = (double[])Best_nodes.y.Clone();
                        }

                        StreamWriter sw2 = new StreamWriter(ss3);

                        for (int ii = 0; ii < No_nonanchor; ii++)
                        {
                            sw2.Write(BEST.x[ii] + " ");
                        }
                        sw2.WriteLine();
                        for (int ii = 0; ii < No_nonanchor; ii++)
                        {
                            sw2.Write(BEST.y[ii] + " ");
                        }
                        sw2.Close();
                        sw = new StreamWriter(ss2);
                        for (int tt = 0; tt < Ntest; tt++)
                        {
                            sw.Write(costs[tt].ToString() + " ");
                        }
                        sw.Close();
                    }
                }
                //}
            }

            #endregion
            #region PMS
            if (method == "PMS")
            {

                double[] h = new double[No_nonanchor];
                for (int ii = 0; ii < No_nonanchor; ii++)
                    h[ii] = 0.1;
                double[,] A = new double[No_nonanchor, No_nonanchor];
                //double[] inheritance_factors = { 0.01, 0.05, 0.1, 0.3, 0.6 };
                //double[] Alphas = {2,30,200 };
                //double[] Bethas = { 0.7,1,5,10,50 };
                bool updated = false;
                double[] rho_s = { 0.01, 0.01, 0.001, 0.001, 0.01, 0.1, 0.7, 0.01, 0.01, 0.01, 0.01, 0.01, 0.001, 0, 0.001, 0.01, 0.1, 0.001, 0.001, 0.01, 0.01, 0.1 };
                double[] inh_rate = { 0.05, 0.01, 0.05, 0.001, 0.001, 0.05, 0.01, 0.001, 0.01, 0.05, 0.3, 0.05, 0.05, 0, 0.3, 0.01, 0.05, 0.001, 0.3, 0.05, 0.1, 0.01 };
                //for (int rh = 0; rh < 5; rh++)
                //    for (int cr = 0; cr < 5; cr++)
                //    {
                double inheritance_factor = 0.9;
                double root = No_nonanchor * inheritance_factor;
                double Cr = 1.0 / (Math.Pow(2, 1.0 / root));
                root = No_nonanchor * (1 - inheritance_factor);
                int local_budget = 150;
                double alpha = 2;
                double betha = 0.5;
                double rho = 0.4;
                string ss2 = "fit_N_Uniform" + No_nonanchor.ToString() + "pop_" + NPop.ToString() + "Localization" + method + "inh_rate_" + inheritance_factor.ToString() + "rho_" + rho.ToString() + "_test_" + Ntest.ToString() + ".txt";
                #region file_existency
                if (!File.Exists(ss2))
                {
                    for (int tt = 0; tt < Ntest; tt++)
                    {
                        rho = 0.4;
                        for (int ii = 0; ii < No_nonanchor; ii++)
                            A[ii, ii] = 1;
                        indi Elite = new indi();
                        indi Indi = new indi();
                        ////////////////////////////////////
                        //////////////////////////////
                        Elite = NFC.Initialization(No_nonanchor, Rand);
                        Elite.cost = NFC.compute_cost(Elite.x, Elite.y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);
                        for (int it = 0; it < Iteration * NPop; )
                        {
                            ///////////////////Long distance/////////////////////////////////
                            updated = false;
                            int global_budget_spent_gen = Convert.ToInt32(0.05 * Iteration * NPop);
                            int gg = 0;
                            while (!updated && it < Iteration * NPop && gg < global_budget_spent_gen)
                            {
                                #region Long_distance
                                int checker = 0;
                                Indi = NFC.Initialization(No_nonanchor, Rand);
                                int i = Rand.Next(No_nonanchor);
                                Indi.x[i] = Elite.x[i]; Indi.y[i] = Elite.y[i];
                                while (Rand.NextDouble() <= Cr)
                                {
                                    Indi.x[i] = Elite.x[i];
                                    Indi.y[i] = Elite.y[i];
                                    i++;
                                    if (i == No_nonanchor)
                                        i = 0;
                                    if (checker >= No_nonanchor)
                                        break;
                                    checker++;
                                }
                                Indi.cost = NFC.compute_cost(Indi.x, Indi.y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);
                                if (Indi.cost <= Elite.cost)
                                {
                                    Elite.x = (double[])Indi.x.Clone();
                                    Elite.y = (double[])Indi.y.Clone();
                                    Elite.cost = Indi.cost;
                                    updated = true;
                                }
                                //////////////////////Function Evaluation increasing///////////////////
                                it++;
                                gg++;
                                #endregion
                            }
                            /////////////////////////Short Distance////////////////////////////////////////
                            double decisional_par = Rand.NextDouble();
                            if (decisional_par < 0.5)
                            {
                                #region Short_region
                                indi indi2 = new indi();
                                for (int ii = 0; ii < local_budget && it < NPop * Iteration; ii++)
                                {
                                    Indi.x = (double[])Elite.x.Clone();
                                    indi2.x = (double[])Elite.x.Clone();
                                    for (int jj = 0; jj < No_nonanchor; jj++)
                                    {
                                        if (Elite.x[jj] - rho > -1)
                                            indi2.x[jj] = Elite.x[jj] - rho;
                                        else
                                            indi2.x[jj] = Elite.x[jj];
                                        indi2.cost = NFC.compute_cost(indi2.x, indi2.y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);
                                        it++;
                                        if (indi2.cost >= Indi.cost)
                                        {
                                            Indi.x = (double[])indi2.x.Clone();
                                            Indi.y = (double[])indi2.y.Clone();
                                            Indi.cost = indi2.cost;
                                        }
                                        else
                                        {
                                            if (Elite.x[jj] + rho / 2 < 1)
                                                indi2.x[jj] = Elite.x[jj] + rho / 2;
                                            else
                                                indi2.x[jj] = Elite.x[jj];
                                            indi2.cost = NFC.compute_cost(indi2.x, indi2.y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);
                                            it++;
                                            if (indi2.cost <= Indi.cost)
                                            {
                                                Indi.x = (double[])indi2.x.Clone();
                                                Indi.y = (double[])indi2.y.Clone();
                                                Indi.cost = indi2.cost;
                                            }
                                        }
                                    }
                                    if (Indi.cost <= Elite.cost)
                                    {
                                        Elite.x = (double[])Indi.x.Clone();
                                        Elite.y = (double[])Indi.y.Clone();
                                        Elite.cost = Indi.cost;
                                        updated = true;
                                    }
                                    else
                                        rho /= 2;

                                }
                                #endregion
                            }
                            else
                            {
                                NFC.rosenbrock_localization(50000, alpha, betha, No_nonanchor, ref Elite.cost, ref Elite.x, ref Elite.y, Rand, range, Realx, Realy, gaussian, Noise_Factor, Anchorx, Anchory, ref it);
                            }
                        }
                        costs[tt] = Elite.cost;
                    }
                    sw = new StreamWriter(ss2);
                    for (int tt = 0; tt < Ntest; tt++)
                    {
                        sw.Write(costs[tt].ToString() + " ");
                    }
                    sw.Close();

                    //}
                }
                #endregion


            }
            #endregion
            #region CCPSO2
            else if (method == "CCPSO2")
            {
                //StreamReader SR = new StreamReader("cauchy.txt");
                //string all = SR.ReadToEnd();
                //string[] all_with_space = all.Split('\n', ',');
                //double[] Cauchy = new double[50000];
                //for (int ii = 0; ii < 50000; ii++)
                //{
                //    Cauchy[ii] = Convert.ToDouble(all_with_space[ii]);
                //}
                double[,] S1_A = { { 0.02, 0.04, 0.08, 0.16, 0.32 }, { 0.02, 0.05, 0.1, 0.5, 1 }, { 0.04, 0.16, 0.32, 0.64, 1 }, { 0.1, 0.3, 0.5, 0.7, 0.9 }, { 0.02, 0.2, 0.4, 0.8, 1 } };
                int[] S = { 2, 5, 10, 50, 100 };
                int[] SSS = { 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 2, 0, 4, 2, 4, 4, 4, 3, 4, 3 };
                double[] p_s = { 0.001, 0.001, 0, 0.1, 0, 0, 0, 0.1, 0.001, 0.001, 0, 0.1, 0.001, 0, 0, 0.001, 0.1, 1, 0, 0.5, 1, 0 };
                //for (int sss = 0; sss < 5; sss++)
                //{
                for (int ii = 0; ii < S.Length; ii++)
                    S[ii] = Convert.ToInt32(S1_A[SSS[4], ii] * No_nonanchor);
                double p = 0;// it has a chance to be a choice as a parameter
                indi[] Pop = new indi[NPop];
                indi BEST = new indi();
                int[] Indexes = new int[No_nonanchor];
                string ss2 = "fit_NUniform_" + No_nonanchor.ToString() + "_Pop_" + NPop.ToString() + "Locali" + "p_" + p.ToString() + "S" + SSS[4].ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "range" + range.ToString() + ".txt";
                string ss3 = "loc_NUniform_" + No_nonanchor.ToString() + "_Pop_" + NPop.ToString() + "Locali" + "p_" + p.ToString() + "S" + SSS[4].ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "range" + range.ToString() + ".txt";
                if (!File.Exists(ss2))
                {
                    for (int tt = 0; tt < Ntest; tt++)
                    {
                        int ind = Rand.Next(S.Length);
                        int s = S[ind];
                        int K = No_nonanchor / s;
                        Swarms_loc[] SWs = new Swarms_loc[K];
                        Swarms_loc[] LocalBest = new Swarms_loc[K];
                        Swarms_loc[] PersonalBest = new Swarms_loc[K];
                        indi[] PersonBest = new indi[NPop];
                        indi Ybest = new indi();
                        Ybest.cost = 1E290;
                        Swarms_loc[] PYbest = new Swarms_loc[K];
                        for (int ii = 0; ii < NPop; ii++)
                        {
                            Pop[ii] = NFC.Initialization(No_nonanchor, Rand);
                            Pop[ii].cost = NFC.compute_cost(Pop[ii].x, Pop[ii].y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);
                            PersonBest[ii].x = (double[])Pop[ii].x.Clone();
                            PersonBest[ii].y = (double[])Pop[ii].y.Clone();
                            if (Ybest.cost > Pop[ii].cost)
                            {
                                Ybest.x = (double[])Pop[ii].x.Clone();
                                Ybest.y = (double[])Pop[ii].y.Clone();
                                Ybest.cost = Pop[ii].cost;
                            }
                            PersonBest[ii].cost = Pop[ii].cost;
                        }
                        double p_yfirst = Ybest.cost;
                        for (int it = 0; it < NPop * Iteration; )
                        {
                            if (p_yfirst >= Ybest.cost && it != 0)
                            {
                                ind = Rand.Next(S.Length);
                                s = S[ind];
                                K = No_nonanchor / s;
                                SWs = new Swarms_loc[K];
                                LocalBest = new Swarms_loc[K];
                                PersonalBest = new Swarms_loc[K];
                                PYbest = new Swarms_loc[K];
                            }
                            p_yfirst = Ybest.cost;
                            Indexes = NFC.sq_rand_gen2(No_nonanchor, No_nonanchor, Rand);
                            int ll = 0;
                            for (int kk = 0; kk < K; kk++)
                            {
                                SWs[kk].particles = new indi[NPop];
                                LocalBest[kk].particles = new indi[NPop];
                                PersonalBest[kk].particles = new indi[NPop];
                                PYbest[kk].particles = new indi[1];
                                for (int jj = 0; jj < NPop; jj++)
                                {
                                    SWs[kk].particles[jj].x = new double[s];
                                    LocalBest[kk].particles[jj].x = new double[s];
                                    PersonalBest[kk].particles[jj].x = new double[s];
                                    PYbest[kk].particles[0].x = new double[s];

                                    SWs[kk].particles[jj].y = new double[s];
                                    LocalBest[kk].particles[jj].y = new double[s];
                                    PersonalBest[kk].particles[jj].y = new double[s];
                                    PYbest[kk].particles[0].y = new double[s];
                                    for (int ii = 0; ii < s; ii++)
                                    {
                                        SWs[kk].particles[jj].x[ii] = Pop[jj].x[Indexes[kk * s + ii]];
                                        PersonalBest[kk].particles[jj].x[ii] = PersonBest[jj].x[Indexes[kk * s + ii]];

                                        SWs[kk].particles[jj].y[ii] = Pop[jj].y[Indexes[kk * s + ii]];
                                        PersonalBest[kk].particles[jj].y[ii] = PersonBest[jj].y[Indexes[kk * s + ii]];
                                    }
                                }
                                for (int ii = 0; ii < s; ii++)
                                {
                                    PYbest[kk].particles[0].x[ii] = Ybest.x[ll];
                                    PYbest[kk].particles[0].y[ii] = Ybest.y[ll++];
                                }

                                PYbest[kk].particles[0].cost = Ybest.cost;
                            }
                            for (int ss = 0; ss < K; ss++)
                            {
                                for (int jj = 0; jj < NPop; jj++)
                                {
                                    double[] constructingx = NFC.Concatenation(Ybest.x, No_nonanchor, ss * s, s, SWs[ss].particles[jj].x);
                                    double[] constructingy = NFC.Concatenation(Ybest.y, No_nonanchor, ss * s, s, SWs[ss].particles[jj].y);
                                    double[] personal_bestx = NFC.Concatenation(Ybest.x, No_nonanchor, ss * s, s, PersonalBest[ss].particles[jj].x);
                                    double[] personal_besty = NFC.Concatenation(Ybest.y, No_nonanchor, ss * s, s, PersonalBest[ss].particles[jj].y);
                                    double[] y_bestx = NFC.Concatenation(Ybest.x, No_nonanchor, ss * s, s, PYbest[ss].particles[0].x);
                                    double[] y_besty = NFC.Concatenation(Ybest.y, No_nonanchor, ss * s, s, PYbest[ss].particles[0].y);
                                    SWs[ss].particles[jj].cost = NFC.compute_cost(constructingx, constructingy, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);
                                    PersonalBest[ss].particles[jj].cost = NFC.compute_cost(personal_bestx, personal_besty, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);
                                    if (PersonalBest[ss].particles[jj].cost > SWs[ss].particles[jj].cost)
                                    {
                                        PersonalBest[ss].particles[jj].cost = SWs[ss].particles[jj].cost;
                                        PersonalBest[ss].particles[jj].x = (double[])SWs[ss].particles[jj].x.Clone();
                                        PersonalBest[ss].particles[jj].y = (double[])SWs[ss].particles[jj].y.Clone();
                                    }
                                    PYbest[ss].particles[0].cost = NFC.compute_cost(y_bestx, y_besty, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);
                                    if (PersonalBest[ss].particles[jj].cost < PYbest[ss].particles[0].cost)
                                    {
                                        PYbest[ss].particles[0].cost = PersonalBest[ss].particles[jj].cost;
                                        PYbest[ss].particles[0].x = (double[])PersonalBest[ss].particles[jj].x.Clone();
                                        PYbest[ss].particles[0].y = (double[])PersonalBest[ss].particles[jj].y.Clone();
                                    }
                                    it = it + 3;
                                }
                                /////////////////////////////////////////////////
                                for (int jj = 0; jj < NPop; jj++)
                                {
                                    int[] Nei = NFC.Ring_withmyself(jj, NPop);
                                    double max = -1E250;
                                    for (int ii = 0; ii < Nei.Length; ii++)
                                    {
                                        if (PersonalBest[ss].particles[Nei[ii]].cost > max)
                                        {
                                            LocalBest[ss].particles[jj].cost = PersonalBest[ss].particles[Nei[ii]].cost;
                                            LocalBest[ss].particles[jj].x = (double[])PersonalBest[ss].particles[Nei[ii]].x.Clone();
                                            LocalBest[ss].particles[jj].y = (double[])PersonalBest[ss].particles[Nei[ii]].y.Clone();
                                            max = PersonalBest[ss].particles[Nei[ii]].cost;
                                        }
                                    }
                                }
                                if (PYbest[ss].particles[0].cost < Ybest.cost)
                                {
                                    Ybest.cost = PYbest[ss].particles[0].cost;
                                    int mm = 0;
                                    for (int ii = ss * s; ii < ss * s + s; ii++)
                                    {
                                        Ybest.x[ii] = PYbest[ss].particles[0].x[mm];
                                        Ybest.y[ii] = PYbest[ss].particles[0].y[mm++];
                                    }
                                }


                            }
                            for (int ss = 0; ss < K; ss++)
                            {
                                for (int jj = 0; jj < NPop; jj++)
                                {
                                    for (int kk = 0; kk < s; kk++)
                                    {
                                        if (Rand.NextDouble() <= p)
                                        {
                                            SWs[ss].particles[jj].x[kk] = PersonalBest[ss].particles[jj].x[kk] + NFC.Cauchy_mu(Rand) * (PersonalBest[ss].particles[jj].x[kk] - LocalBest[ss].particles[jj].x[kk]);
                                            SWs[ss].particles[jj].y[kk] = PersonalBest[ss].particles[jj].y[kk] + NFC.Cauchy_mu(Rand) * (PersonalBest[ss].particles[jj].y[kk] - LocalBest[ss].particles[jj].y[kk]);
                                        }
                                        else
                                        {
                                            SWs[ss].particles[jj].x[kk] = LocalBest[ss].particles[jj].x[kk] + NFC.Normal_Distribution(Rand, 0, 1) * (PersonalBest[ss].particles[jj].x[kk] - LocalBest[ss].particles[jj].x[kk]);
                                            SWs[ss].particles[jj].y[kk] = LocalBest[ss].particles[jj].y[kk] + NFC.Normal_Distribution(Rand, 0, 1) * (PersonalBest[ss].particles[jj].y[kk] - LocalBest[ss].particles[jj].y[kk]);
                                        }
                                        Pop[jj].x[ss * s + kk] = SWs[ss].particles[jj].x[kk];
                                        PersonBest[jj].x[ss * s + kk] = PersonalBest[ss].particles[jj].x[kk];
                                        Pop[jj].y[ss * s + kk] = SWs[ss].particles[jj].y[kk];
                                        PersonBest[jj].y[ss * s + kk] = PersonalBest[ss].particles[jj].y[kk];
                                    }
                                }
                            }
                            for (int jj = 0; jj < NPop; jj++)
                            {
                                NFC.Node_treatment_localization(ref Pop[jj].x, ref Pop[jj].y, No_nonanchor, Rand);
                            }

                        }

                        costs[tt] = NFC.Calc_error(Ybest.x, Ybest.y, Realx, Realy, range);
                        BEST.x = (double[])Ybest.x.Clone();
                        BEST.y = (double[])Ybest.y.Clone();
                    }
                    StreamWriter sw2 = new StreamWriter(ss3);

                    for (int ii = 0; ii < No_nonanchor; ii++)
                    {
                        sw2.Write(BEST.x[ii] + " ");
                    }
                    sw2.WriteLine();
                    for (int ii = 0; ii < No_nonanchor; ii++)
                    {
                        sw2.Write(BEST.y[ii] + " ");
                    }
                    sw2.Close();
                    sw = new StreamWriter(ss2);
                    for (int tt = 0; tt < Ntest; tt++)
                    {
                        sw.Write(costs[tt].ToString() + " ");
                    }
                    sw.Close();
                }
                //}
                //}
            }

            #endregion
            #region JADE
            else if (method == "JADE")
            {

                double[] c_a = { 0.1, 0.1, 0.0001, 0, 0.0001, 0, 0, 0.1, 0, 0.0001, 0, 0.0001, 0.0001, 0, 0, 0.001, 0, 0.001, 0.0001, 0.01, 0, 0.001 };
                double[] p_a = { 0.01, 0.3, 0.3, 0.3, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.3, 0.05, 0.05, 0, 0.05, 0.5, 0.05, 0.1, 0.5, 0.5, 0.05, 0.1 };//greediness

                //double[] Fs = { 0.1, 0.5, 0.5, 0.1, 0.1, 0.1, 0.5, 0.5, 0.1, 0.5, 0.1, 0.1, 0.1,0, 0.1, 0.1, 0.1, 0.5, 2.0, 0.1, 0.1, 0.1 };
                //double[] C_s = { 0.2, 0.2, 0.2, 0.4, 0.2, 0.4, 0.2, 0.2, 0.4, 0.2, 0.2, 0.2, 0.2, 0,0.4, 1, 1, 0.2, 0.2, 1, 0.4, 0.8 };
                double[] Cr_A = new double[NPop];
                double[] F_A = new double[NPop];
                double Cr = 0.5;
                double F = 0.5;
                //for (int CC = 0; CC < 5; CC++)
                //    for (int PP = 0; PP < 5; PP++)
                //    {
                double c = 0.0001;
                double p = 0.1;
                int p_group = Convert.ToInt32(100 * p % NPop);
                //CD.Scale=0.1;
                indi BEST = new indi();
                string ss2 = "fit_N_Uniform" + No_nonanchor.ToString() + "_Pop_" + NPop.ToString() + "_c_" + c.ToString() + "_p_" + p.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "Range" + range.ToString() + ".txt";
                string ss3 = "loc_N_Uniform" + No_nonanchor.ToString() + "_Pop_" + NPop.ToString() + "_c_" + c.ToString() + "_p_" + p.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "Range" + range.ToString() + ".txt";
                if (!File.Exists(ss2))
                {
                    for (int tt = 0; tt < Ntest; tt++)
                    {
                        Cr_A = new double[NPop];
                        F_A = new double[NPop];
                        Cr = 0.5;
                        F = 0.5;
                        indi[] Pop = new indi[NPop];
                        indi[] V = new indi[NPop];
                        indi[] U = new indi[NPop];
                        indi[] Best = new indi[NPop];
                        ArrayList SF = new ArrayList();
                        ArrayList Scr = new ArrayList();
                        ArrayList A_x = new ArrayList();
                        ArrayList A_y = new ArrayList();
                        ArrayList A_f = new ArrayList();
                        indi Best_nodes = new indi();
                        double[] distances = new double[NPop];
                        Best_nodes.cost = 1E250;
                        for (int ii = 0; ii < NPop; ii++)
                        {
                            Pop[ii] = NFC.Initialization(No_nonanchor, Rand);
                            Pop[ii].cost = NFC.compute_cost(Pop[ii].x, Pop[ii].y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);
                            A_x.Add(Pop[ii].x.Clone());
                            A_y.Add(Pop[ii].y.Clone());
                            A_f.Add(Pop[ii].cost);
                            /////////init V
                            V[ii].x = new double[No_nonanchor];
                            U[ii].x = new double[No_nonanchor];
                            V[ii].y = new double[No_nonanchor];
                            U[ii].y = new double[No_nonanchor];
                        }
                        ///////////////////////////////////////////////////////////////////////////////////////
                        for (int it = 0; it < Iteration; it++)
                        {
                            Best = NFC.mergesort(Pop, NPop);
                            for (int ii = 0; ii < NPop; ii++)
                            {
                                F_A[ii] = NFC.Cauchy(Rand, F, 0.1);
                                Cr_A[ii] = NFC.Normal_Distribution(Rand, Cr, 0.1);
                                int best_index = Rand.Next(p_group);
                                indi X_best = Best[best_index];
                                int index_r1 = Rand.Next(NPop);
                                if (X_best.cost == Pop[index_r1].cost)
                                {
                                    while (NFC.similarity(X_best.x, Pop[index_r1].x) && NFC.similarity(X_best.y, Pop[index_r1].y))
                                        index_r1 = Rand.Next(NPop);
                                }
                                int index_r2 = Rand.Next(A_f.Count);
                                double Fit_A = Convert.ToDouble(A_f[index_r2]);
                                if (X_best.cost == Fit_A || Fit_A == Pop[index_r1].cost)
                                {
                                    double[] X = (double[])A_x[index_r2];
                                    double[] Y = (double[])A_y[index_r2];
                                    while (NFC.similarity(X_best.x, X) && NFC.similarity(X_best.y, Y))
                                    {
                                        index_r2 = Rand.Next(NPop);
                                        X = (double[])A_x[index_r2];
                                        Y = (double[])A_y[index_r2];
                                    }
                                }
                                int J_rand = Rand.Next(No_nonanchor);
                                for (int jj = 0; jj < No_nonanchor; jj++)
                                {
                                    double[] x = (double[])A_x[index_r2];
                                    V[ii].x[jj] = Pop[ii].x[jj] + F_A[ii] * (X_best.x[jj] - Pop[ii].x[jj]) + F_A[ii] * (Pop[index_r1].x[jj] - x[jj]);
                                    double[] y = (double[])A_y[index_r2];
                                    V[ii].y[jj] = Pop[ii].y[jj] + F_A[ii] * (X_best.y[jj] - Pop[ii].y[jj]) + F_A[ii] * (Pop[index_r1].y[jj] - y[jj]);
                                    if ((jj == J_rand) || (Rand.NextDouble() < Cr_A[ii]))
                                    {
                                        U[ii].x[jj] = V[ii].x[jj];
                                        U[ii].y[jj] = V[ii].y[jj];

                                    }
                                    else
                                    {
                                        U[ii].x[jj] = Pop[ii].x[jj];
                                        U[ii].y[jj] = Pop[ii].y[jj];
                                    }
                                }
                                NFC.Node_treatment_localization(ref U[ii].x, ref U[ii].y, No_nonanchor, Rand);
                                U[ii].cost = NFC.compute_cost(U[ii].x, U[ii].y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);
                                if (U[ii].cost <= Pop[ii].cost)
                                {
                                    int rand = Rand.Next(NPop);
                                    A_x.Add((double[])Pop[ii].x.Clone());
                                    A_y.Add((double[])Pop[ii].y.Clone());
                                    A_f.Add(Pop[ii].cost);
                                    Pop[ii].x = (double[])U[ii].x.Clone();
                                    Pop[ii].y = (double[])U[ii].y.Clone();
                                    Pop[ii].cost = U[ii].cost;
                                    SF.Add(F_A[ii]);
                                    Scr.Add(Cr_A[ii]);
                                }
                                if (Best_nodes.cost > Pop[ii].cost)
                                {
                                    Best_nodes.cost = Pop[ii].cost;
                                    Best_nodes.x = (double[])Pop[ii].x.Clone();
                                    Best_nodes.y = (double[])Pop[ii].y.Clone();

                                }
                            }

                            int Overhead = A_x.Count - NPop;
                            for (int ii = 0; ii < Overhead; ii++)
                            {
                                int rand = Rand.Next(A_x.Count);
                                A_x.RemoveAt(rand);
                                A_y.RemoveAt(rand);
                                A_f.RemoveAt(rand);
                            }

                            Cr = (1 - c) * Cr + c * NFC.mean_simple(Scr);
                            F = (1 - c) * F + c * NFC.mean_Lehmer(SF);
                        }

                        costs[tt] = NFC.Calc_error(Best_nodes.x, Best_nodes.y, Realx, Realy, range);
                        BEST.x = (double[])Best_nodes.x.Clone();
                        BEST.y = (double[])Best_nodes.y.Clone();
                    }
                    StreamWriter sw2 = new StreamWriter(ss3);

                    for (int ii = 0; ii < No_nonanchor; ii++)
                    {
                        sw2.Write(BEST.x[ii] + " ");
                    }
                    sw2.WriteLine();
                    for (int ii = 0; ii < No_nonanchor; ii++)
                    {
                        sw2.Write(BEST.y[ii] + " ");
                    }
                    sw2.Close();
                    sw = new StreamWriter(ss2);
                    for (int tt = 0; tt < Ntest; tt++)
                    {
                        sw.Write(costs[tt].ToString() + " ");
                    }
                    sw.Close();
                }
                //}
            }

            #endregion
            #region FES
            else if (method == "FES")
            {

                //double[] Lambda_s = { 0.2, 0.4, 0.6, 0.8, 1 };
                //// sigma constant which is the constant factor to be multiplied to sigma
                //double[] SC_s = { 1.1, 1.2, 1.5, 2, 5 };
                //for (int cc = 0; cc < 5; cc++)
                //    for (int ww = 0; ww < 5; ww++)
                //    {
                indi BEST = new indi();
                double[] Lambda_s = { 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 0, 1.0000, 0.8000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 };
                double[] SC_s = { 1.1000, 1.1000, 1.1000, 1.1000, 1.1000, 1.1000, 1.1000, 1.1000, 1.1000, 1.1000, 1.1000, 1.1000, 1.1000, 0, 1.1000, 2.0000, 1.1000, 1.1000, 1.2000, 1.1000, 1.1000, 1.1000 };
                int breeds_number = 1;
                double SC = 1.1;
                //ss = "Result_N_" + N.ToString() + "_Pop_" + NPop.ToString() + "_Lambda_" + Lambda_s[cc].ToString() + "_Sigma_Fac_" + SC.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "_OF" + OF.ToString() + ".txt";
                string ss2 = "fit_N_Uniform" + No_nonanchor.ToString() + "_Pop_" + NPop.ToString() + "_Lambda_" + breeds_number.ToString() + "_Sigma_Fac_" + SC.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "Range" + range.ToString() + ".txt";
                string ss3 = "fit_N_Uniform" + No_nonanchor.ToString() + "_Pop_" + NPop.ToString() + "_Lambda_" + breeds_number.ToString() + "_Sigma_Fac_" + SC.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "Range" + range.ToString() + ".txt";
                if (!File.Exists(ss2))
                {
                    for (int tt = 0; tt < Ntest; tt++)
                    {
                        indi[] Pop = new indi[NPop];
                        indi Best_nodes = new indi();
                        Best_nodes.cost = 1E250;
                        indi[] Breeds = new indi[breeds_number];
                        double[] distances = new double[NPop];
                        int phi_count = 0;
                        double sigma = 1;
                        for (int ii = 0; ii < NPop; ii++)
                        {
                            Pop[ii] = NFC.Initialization(No_nonanchor, Rand);
                            Pop[ii].cost = NFC.compute_cost(Pop[ii].x, Pop[ii].y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);
                        }

                        for (int it = 0; it < Iteration; it++)
                        {
                            indi[] All = new indi[NPop + breeds_number];
                            for (int ii = 0; ii < breeds_number; ii++)
                            {
                                Breeds[ii].x = new double[No_nonanchor];
                                Breeds[ii].y = new double[No_nonanchor];
                                int parent1 = Rand.Next(NPop);
                                int parent2 = Rand.Next(NPop);
                                for (int n = 0; n < No_nonanchor; n++)
                                {
                                    double rand = Rand.NextDouble();
                                    if (rand < 0.5)
                                    {
                                        Breeds[ii].x[n] = Pop[parent1].x[n];
                                        Breeds[ii].y[n] = Pop[parent1].y[n];
                                    }
                                    else
                                    {
                                        Breeds[ii].x[n] = Pop[parent2].x[n];
                                        Breeds[ii].y[n] = Pop[parent2].y[n];
                                    }
                                }
                            }
                            for (int ii = 0; ii < breeds_number; ii++)
                            {
                                for (int n = 0; n < No_nonanchor; n++)
                                {
                                    Breeds[ii].x[n] = Breeds[ii].x[n] + sigma * NFC.Cauchy(Rand, 0, 1);
                                    Breeds[ii].y[n] = Breeds[ii].y[n] + sigma * NFC.Cauchy(Rand, 0, 1);
                                }
                                NFC.Node_treatment_localization(ref Breeds[ii].x, ref Breeds[ii].y, No_nonanchor, Rand);
                                Breeds[ii].cost = NFC.compute_cost(Breeds[ii].x, Breeds[ii].y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);
                                All[ii].x = (double[])Breeds[ii].x.Clone();
                                All[ii].y = (double[])Breeds[ii].y.Clone();
                                All[ii].cost = Breeds[ii].cost;

                            }

                            double[] OldCost = new double[NPop];
                            for (int ii = 0; ii < NPop; ii++)
                            {
                                OldCost[ii] = NFC.compute_cost(Pop[ii].x, Pop[ii].y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);

                                for (int n = 0; n < No_nonanchor; n++)
                                {
                                    Pop[ii].x[n] = Pop[ii].x[n] + sigma * NFC.Cauchy(Rand, 0, 1);
                                    Pop[ii].y[n] = Pop[ii].y[n] + sigma * NFC.Cauchy(Rand, 0, 1);
                                }
                                NFC.Node_treatment_localization(ref Pop[ii].x, ref Pop[ii].y, No_nonanchor, Rand);
                                Pop[ii].cost = NFC.compute_cost(Pop[ii].x, Pop[ii].y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);
                                // putting all pop into a temporary combinator
                                All[ii + breeds_number].x = (double[])Pop[ii].x.Clone();
                                All[ii + breeds_number].y = (double[])Pop[ii].y.Clone();
                                All[ii + breeds_number].cost = Pop[ii].cost;

                                if (Pop[ii].cost < OldCost[ii])
                                    phi_count++;
                            }
                            int mutation_count = it * NPop;
                            if (phi_count > (mutation_count / 5))
                                sigma = sigma / SC;
                            else if (phi_count < (mutation_count / 5))
                                sigma = sigma * SC;
                            // Sort all POP+Breed
                            All = NFC.mergesort(All, All.Length);

                            //storing the best node into best node container
                            if (Best_nodes.cost > All[0].cost)
                            {
                                Best_nodes.x = (double[])All[0].x.Clone();
                                Best_nodes.y = (double[])All[0].y.Clone();
                                Best_nodes.cost = All[0].cost;
                            }

                            for (int ii = 0; ii < NPop; ii++)
                            {
                                Pop[ii].x = (double[])All[ii].x.Clone();
                                Pop[ii].y = (double[])All[ii].y.Clone();
                                Pop[ii].cost = All[ii].cost;
                            }
                        }

                        costs[tt] = NFC.Calc_error(Best_nodes.x, Best_nodes.y, Realx, Realy, range);
                        BEST.x = (double[])Best_nodes.x.Clone();
                        BEST.y = (double[])Best_nodes.y.Clone();
                    }
                    sw = new StreamWriter(ss3);
                    for (int tt = 0; tt < No_nonanchor; tt++)
                    {
                        sw.Write(BEST.x[tt].ToString() + " ");
                    }
                    sw.WriteLine();
                    for (int tt = 0; tt < No_nonanchor; tt++)
                    {
                        sw.Write(BEST.y[tt].ToString() + " ");
                    }
                    sw.Close();
                    sw = new StreamWriter(ss2);
                    for (int tt = 0; tt < Ntest; tt++)
                    {
                        sw.Write(costs[tt].ToString() + " ");
                    }
                    sw.Close();
                }
                //}

            }

            #endregion
            #region PSO
            else if (method == "PSO")
            {
                double W1 = Convert.ToDouble(Wmin.Text);
                double W2 = Convert.ToDouble(Wmax.Text);
                double w_step = Convert.ToDouble(W_step.Text);
                double C11 = Convert.ToDouble(cmin.Text);
                double C12 = Convert.ToDouble(cmax.Text);
                double C11_step = Convert.ToDouble(c_step.Text);
                double C21 = Convert.ToDouble(cmin21.Text);
                double C22 = Convert.ToDouble(cmin22.Text);
                double C22_step = Convert.ToDouble(c2_step.Text);
                double[] Cs = { 0.0001, 0.001, 0.1, 0.5, 1, 1.5 };
                double[] Ws = { 0.1, 0.5, 0.73, 1, 1.5 };
                //for (double C1 = C11; C1 <= C12; C1 += C11_step)
                //    for (double C2 = C21; C2 <= C12; C2 += C22_step)
                //for (int WW = 0; WW < 5; WW++)
                //    for (int CC = 0; CC < 6; CC++)
                //    {
                //double[] C_s = { 0.5, 0.001, 0.001, 0.001, 0.001, 0.001, 1, 1, 0.001, 0.001, 0.001, 0.001, 0.001, 0, 0.001, 1, 1, 1, 1, 1, 0.001, 1 };
                //double[] W_s = { 0.5, 1, 1, 1, 1, 1, 0.5, 0.5, 1, 1, 1, 1, 1, 0, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1 };
                double C1 = 0.01;
                double C2 = 0.01;
                double W = 1;
                indi BEST = new indi();
                string ss2 = "fit_N_Uniform" + No_nonanchor.ToString() + "_Pop_" + NPop.ToString() + "_W_" + W.ToString() + "_C1_" + C1.ToString() + "_C2_" + C2.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "range" + range.ToString() + ".txt";
                string ss3 = "loc_N_Uniform" + No_nonanchor.ToString() + "_Pop_" + NPop.ToString() + "_W_" + W.ToString() + "_C1_" + C1.ToString() + "_C2_" + C2.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "range" + range.ToString() + ".txt";
                if (!File.Exists(ss2))
                {
                    for (int tt = 0; tt < Ntest; tt++)
                    {
                        indi[] Nodes = new indi[NPop];
                        indi[] BNodes = new indi[NPop];
                        indi[] Forces = new indi[NPop];
                        indi[] a = new indi[NPop];
                        indi[] v = new indi[NPop];
                        indi Best_nodes = new indi();
                        double[] distances = new double[NPop];
                        Best_nodes.cost = 1E290;

                        for (int ii = 0; ii < NPop; ii++)
                        {
                            v[ii].x = new double[No_nonanchor];
                            v[ii].y = new double[No_nonanchor];
                        }
                        for (int ii = 0; ii < NPop; ii++)
                        {
                            Nodes[ii] = NFC.Initialization(No_nonanchor, Rand);
                            Nodes[ii].cost = NFC.compute_cost(Nodes[ii].x, Nodes[ii].y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);
                            BNodes[ii].x = (double[])Nodes[ii].x.Clone();
                            BNodes[ii].y = (double[])Nodes[ii].y.Clone();
                            BNodes[ii].cost = Nodes[ii].cost;
                        }
                        for (int it = 0; it < Iteration; it++)
                        {
                            double Fmin = 1E290; int min_index = 0;
                            for (int ii = 0; ii < NPop; ii++)
                            {
                                Nodes[ii].cost = NFC.compute_cost(Nodes[ii].x, Nodes[ii].y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);
                                if (Fmin > Nodes[ii].cost)
                                {
                                    Fmin = Nodes[ii].cost;
                                    min_index = ii;
                                }
                            }
                            if (Best_nodes.cost > Fmin)
                            {
                                Best_nodes.cost = Fmin;
                                Best_nodes.x = (double[])Nodes[min_index].x.Clone();
                                Best_nodes.y = (double[])Nodes[min_index].y.Clone();

                            }
                            BNodes = NFC.define_pbest(Nodes, BNodes, NPop);
                            for (int ii = 0; ii < NPop; ii++)
                            {
                                for (int n = 0; n < No_nonanchor; n++)
                                {
                                    v[ii].x[n] = W * v[ii].x[n] + C1 * Rand.NextDouble() * (BNodes[ii].x[n] - Nodes[ii].x[n]) + C2 * Rand.NextDouble() * (Best_nodes.x[n] - Nodes[ii].x[n]);
                                    v[ii].y[n] = W * v[ii].y[n] + C1 * Rand.NextDouble() * (BNodes[ii].y[n] - Nodes[ii].y[n]) + C2 * Rand.NextDouble() * (Best_nodes.y[n] - Nodes[ii].y[n]);
                                    Nodes[ii].x[n] = Nodes[ii].x[n] + v[ii].x[n];
                                    Nodes[ii].y[n] = Nodes[ii].y[n] + v[ii].y[n];
                                }

                                NFC.Node_treatment_localization(ref Nodes[ii].x, ref Nodes[ii].y, No_nonanchor, Rand);
                                NFC.Node_treatment_localization(ref v[ii].x, ref v[ii].y, No_nonanchor, Rand);
                            }
                        }

                        costs[tt] = NFC.Calc_error(Best_nodes.x, Best_nodes.y, Realx, Realy, range);
                        BEST.x = (double[])Best_nodes.x.Clone();
                        BEST.y = (double[])Best_nodes.y.Clone();
                    }


                    sw = new StreamWriter(ss2);
                    StreamWriter sw2 = new StreamWriter(ss3);

                    for (int ii = 0; ii < No_nonanchor; ii++)
                    {
                        sw2.Write(BEST.x[ii] + " ");
                    }
                    sw2.WriteLine();
                    for (int ii = 0; ii < No_nonanchor; ii++)
                    {
                        sw2.Write(BEST.y[ii] + " ");
                    }
                    sw2.Close();
                    for (int tt = 0; tt < Ntest; tt++)
                    {
                        sw.Write(costs[tt].ToString() + " ");
                    }
                    sw.Close();
                }
                //}
            }

            #endregion
            #region GA
            else if (method == "Genetic")
            {
                //    double[] mus = { 0.002, 0.003, 0.004, 0.005, 0.01 };
                double[] Crs = { 0.2, 0.4, 0.6, 0.8, 1 };
                //    for (int muindex = 0; muindex < mus.Length; muindex++)
                //for (int cr_index = 0; cr_index < Crs.Length; cr_index++)
                //{
                //double[] mus = { 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.005, 0.003, 0.003, 0.003, 0.005, 0.003, 0.003,0, 0.003, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.005 };
                //double[] Crs = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.8, 0.4, 0.4, 1.0, 0.8, 0.2, 1.0 };
                double cr_rate = 1;
                double mu_rate = 0.03;
                //double mu_rate = muindex[Convert.ToInt32(MM)];

                indi BEST = new indi();
                string ss2 = "fit_N_Uniform" + No_nonanchor.ToString() + "pop_" + NPop.ToString() + method + "_cr_" + cr_rate.ToString() + "_mu_rate_" + mu_rate.ToString() + "_test_" + Ntest.ToString() + "range" + range.ToString() + ".txt";
                string ss3 = "local_Uniform_" + No_nonanchor.ToString() + "pop_" + NPop.ToString() + method + "_cr_" + cr_rate.ToString() + "_mu_rate_" + mu_rate.ToString() + "_test_" + Ntest.ToString() + "range" + range.ToString() + ".txt";
                //string ss2 = "fit_N_" + No_nonanchor.ToString() + "pop_" + NPop.ToString() + method + "_cr_" + cr_rate.ToString() + "_mu_rate_" + mu_rate.ToString() + "_test_" + Ntest.ToString() + "range" + range.ToString() + ".txt";
                //string ss3 = "local_" + No_nonanchor.ToString() + "pop_" + NPop.ToString() + method + "_cr_" + cr_rate.ToString() + "_mu_rate_" + mu_rate.ToString() + "_test_" + Ntest.ToString() + "range" + range.ToString() + ".txt";   
                #region file_existency
                if (!File.Exists(ss2))
                {
                    for (int tt = 0; tt < Ntest; tt++)
                    {
                        indi[] GA = new indi[NPop];
                        double Fmax = Double.NegativeInfinity; double Fmin = 100000000;
                        indi Best_indi = new indi();
                        Best_indi.cost = Double.PositiveInfinity;
                        for (int ii = 0; ii < NPop; ii++)
                        {

                            GA[ii] = NFC.Initialization(No_nonanchor, Rand);
                            GA[ii].cost = NFC.compute_cost(GA[ii].x, GA[ii].y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);
                        }
                        for (int it = 0; it < Iteration; it++)
                        {
                            Fmax = Double.NegativeInfinity; Fmin = 100000000;
                            int min_index = -1;
                            for (int ii = 0; ii < NPop; ii++)
                            {
                                GA[ii].cost = NFC.compute_cost(GA[ii].x, GA[ii].y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);
                                if (Fmax < GA[ii].cost)
                                {
                                    Fmax = GA[ii].cost;

                                }
                                if (Fmin > GA[ii].cost)
                                {
                                    Fmin = GA[ii].cost;
                                    min_index = ii;
                                }
                            }
                            if (Fmin == Fmax)
                            {
                                for (int ii = 0; ii < NPop; ii++)
                                {
                                    GA[ii] = NFC.Initialization(No_nonanchor, Rand);
                                    GA[ii].cost = NFC.compute_cost(GA[ii].x, GA[ii].y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);
                                    if (Fmax < GA[ii].cost)
                                        Fmax = GA[ii].cost;
                                    if (Fmin > GA[ii].cost)
                                    {
                                        Fmin = GA[ii].cost;
                                        min_index = ii;
                                    }
                                }
                            }
                            if (Fmin < Best_indi.cost)
                            {
                                Best_indi.cost = Fmin;
                                Best_indi.x = (double[])GA[min_index].x.Clone();
                                Best_indi.y = (double[])GA[min_index].y.Clone();
                            }
                            double[] help = new double[NPop];
                            for (int ii = 0; ii < NPop; ii++)
                            {
                                help[ii] = ((Fmax + Fmin) - GA[ii].cost);
                            }
                            for (int ii = 0; ii < NPop; ii++)
                            {
                                if (cr_rate > Rand.NextDouble())
                                {
                                    int index = NFC.Select_with_Probability(help, Fmin, NPop, Rand);
                                    int index2 = NFC.Select_with_Probability(help, Fmin, NPop, Rand);
                                    indi childs = NFC.Cross_Over_one_point(GA[index], GA[index2], No_nonanchor, Rand, Anchorx, Anchory, range, Realx, Realy, Noise_Factor, gaussian);
                                    GA[ii].x = (double[])childs.x.Clone();
                                    GA[ii].y = (double[])childs.y.Clone();
                                }
                                NFC.mutation_location(ref GA[ii], mu_rate, No_nonanchor, Rand);
                            }
                        }
                        costs[tt] = NFC.Calc_error(Best_indi.x, Best_indi.y, Realx, Realy, range);
                        BEST.x = (double[])Best_indi.x.Clone();
                        BEST.y = (double[])Best_indi.y.Clone();
                    }


                    sw = new StreamWriter(ss2);
                    StreamWriter sw2 = new StreamWriter(ss3);

                    for (int ii = 0; ii < No_nonanchor; ii++)
                    {
                        sw2.Write(Realx[ii] + " ");
                    }
                    sw2.WriteLine();
                    for (int ii = 0; ii < No_nonanchor; ii++)
                    {
                        sw2.Write(Realy[ii] + " ");
                    }
                    sw2.Close();
                    for (int tt = 0; tt < Ntest; tt++)
                    {
                        sw.Write(costs[tt].ToString() + " ");
                    }
                    sw.Close();
                    //}
                }
                #endregion



            }
            #endregion
            #region OMOA
            if (method == "OMOA")
            {
                double Jr = Convert.ToDouble(Jumping_rate.Text);
                double intensity = 1;
                double alpha = 0.6;
                Jr = 0.3;
                string structure = "Cellular";
                string ss = "Result_Pop_" + NPop.ToString() + ".txt";
                if (!File.Exists(ss))
                {
                    for (int tt = 0; tt < Ntest; tt++)
                    {
                        indi[] Nodes = new indi[NPop];
                        indi[] ONodes = new indi[NPop];
                        indi[] Forces = new indi[NPop];
                        indi[] a = new indi[NPop];
                        indi[] v = new indi[NPop];
                        double[] Magnet = new double[NPop];
                        double[] Mass = new double[NPop];
                        indi Best_indi = new indi();
                        double[] distances = new double[NPop];
                        Best_indi.cost = -1E290;
                        for (int ii = 0; ii < NPop; ii++)
                        {
                            Forces[ii].x = new double[No_nonanchor];
                            a[ii].x = new double[No_nonanchor];
                            v[ii].x = new double[No_nonanchor];
                            Forces[ii].y = new double[No_nonanchor];
                            a[ii].y = new double[No_nonanchor];
                            v[ii].y = new double[No_nonanchor];
                        }
                        for (int ii = 0; ii < NPop; ii++)
                        {
                            Nodes[ii] = NFC.Initialization(No_nonanchor, Rand);
                            ONodes[ii].x = NFC.Opposition_Learning(Nodes[ii].x);
                            ONodes[ii].y = NFC.Opposition_Learning(Nodes[ii].y);
                        }
                        for (int it = 0; it < Iteration; it++)
                        {
                            double Fmax = -1E290; double Fmin = 1000000000;
                            int min_index = 0; int max_index = 0;
                            for (int ii = 0; ii < NPop; ii++)
                            {
                                Nodes[ii].cost = NFC.compute_cost(Nodes[ii].x, Nodes[ii].y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);

                                if (Jr > Rand.NextDouble())
                                {

                                    ONodes[ii].x = NFC.generating_jumping(Nodes, NPop, No_nonanchor, ii);
                                    ONodes[ii].y = NFC.generating_jumping(Nodes, NPop, No_nonanchor, ii);
                                    ONodes[ii].cost = NFC.compute_cost(Nodes[ii].x, Nodes[ii].y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);
                                    if (Nodes[ii].cost < ONodes[ii].cost)
                                    {
                                        Nodes[ii].x = (double[])ONodes[ii].x.Clone();
                                        Nodes[ii].y = (double[])ONodes[ii].y.Clone();
                                        Nodes[ii].cost = ONodes[ii].cost;
                                    }
                                    NFC.Node_treatment_localization(ref Nodes[ii].x, ref Nodes[ii].y, No_nonanchor, Rand);
                                    NFC.Node_treatment_localization(ref v[ii].x, ref v[ii].y, No_nonanchor, Rand);
                                }

                                if (Fmax < Nodes[ii].cost)
                                {
                                    Fmax = Nodes[ii].cost;
                                    max_index = ii;
                                }
                                if (Fmin > Nodes[ii].cost)
                                {
                                    Fmin = Nodes[ii].cost;
                                    min_index = ii;
                                }
                            }

                            if (Best_indi.cost < Fmax)
                            {
                                Best_indi.cost = Fmax;
                                Best_indi.x = (double[])Nodes[max_index].x.Clone();
                                Best_indi.y = (double[])Nodes[max_index].x.Clone();
                            }
                            double ranges = Fmax - Fmin;
                            if (ranges == 0)
                            {
                                Fmax = -1000000000; Fmin = 1000000000;
                                min_index = 0; max_index = 0;
                                for (int ii = 0; ii < NPop; ii++)
                                {
                                    Nodes[ii] = NFC.Initialization(No_nonanchor, Rand);
                                    ONodes[ii].x = NFC.Opposition_Learning(Nodes[ii].x);
                                    ONodes[ii].y = NFC.Opposition_Learning(Nodes[ii].y);
                                }
                            }
                            for (int ii = 0; ii < NPop; ii++)
                            {
                                Magnet[ii] = (Nodes[ii].cost - Fmin) / (range);
                                Mass[ii] = Magnet[ii] * intensity + alpha;
                            }
                            for (int ii = 0; ii < NPop; ii++)
                            {
                                Forces[ii].x = new double[No_nonanchor];
                                Forces[ii].y = new double[No_nonanchor];
                            }
                            for (int ii = 0; ii < NPop; ii++)
                            {
                                int[] Neigbour = NFC.Structures(ii, NPop, structure);
                                for (int jj = 0; jj < Neigbour.Length; jj++)
                                {
                                    if (ii != Neigbour[jj])
                                    {

                                        for (int n = 0; n < No_nonanchor; n++)
                                        {
                                            Forces[ii].x[n] = Forces[ii].x[n] + (Nodes[Neigbour[jj]].x[n] - Nodes[ii].x[n]) * Magnet[Neigbour[jj]] / distances[ii];
                                            Forces[ii].y[n] = Forces[ii].y[n] + (Nodes[Neigbour[jj]].y[n] - Nodes[ii].y[n]) * Magnet[Neigbour[jj]] / distances[ii];
                                        }



                                    }
                                }
                            }
                            for (int ii = 0; ii < NPop; ii++)
                            {
                                for (int n = 0; n < No_nonanchor; n++)
                                {
                                    v[ii].x[n] = (Forces[ii].x[n] / Mass[ii]);
                                    v[ii].y[n] = (Forces[ii].y[n] / Mass[ii]);
                                    Nodes[ii].x[n] = Nodes[ii].x[n] + v[ii].x[n];
                                    Nodes[ii].y[n] = Nodes[ii].y[n] + v[ii].y[n];
                                }

                                NFC.Node_treatment_localization(ref Nodes[ii].x, ref Nodes[ii].y, No_nonanchor, Rand);
                                NFC.Node_treatment_localization(ref v[ii].x, ref v[ii].y, No_nonanchor, Rand);
                            }

                        }
                        costs[tt] = Best_indi.cost;
                    }
                    sw = new StreamWriter(ss);
                    for (int tt = 0; tt < Ntest; tt++)
                    {
                        sw.Write(costs[tt].ToString() + " ");
                    }
                    sw.Close();
                }
                //}
                //}
            }
            #endregion
            #region MASSW_structure
            else if (method == "MASSW_S")
            {
                int O = NPop * Convert.ToInt32(NI.Text);
                //double[] N_LS_GL = { 0.1, 0.3, 0.5, 0.7, 0.9 };
                //double[] Mu_rate_s = {0.00001,0.001,0.01,0.1,0.5 };
                double[] N_LS_GL = { 0.5, 0.7 };
                double[] Mu_rate_s = { 0.01, 0.1, 0.5 };
                string[] Structs = { "Ring", "Star", "Btree", "Bipartie", "Cluster", "Grid", "Ladder", "Cross-Ladder", "Cellular", "basic" };
                //string[] Structs = { "Cellular" };
                //for (int Portion_index = 0; Portion_index < 5; Portion_index++)
                //    for (int M_index = 0; M_index < 5; M_index++)
                //    {
                for (int SS = 0; SS < Structs.Length; SS++)
                {
                    string Struct = Structs[SS];
                    //double[] N_LS_GL = { 0.3,0.5,0.5,0.3,0.5,0.5,0.5,0.5,0.3,0.5,0.7,0.3,0.5,0,0.3,0.5,0.5,0.5,0.5,0.7,0.3,0.5 };
                    //double[] Mu_rate_s = { 5, 0.5, 0.001, 0.001, 5, 5, 0.5, 0.5, 0.001, 0.01, 0.001, 0.01, 0.5, 0, 0.001, 0.1, 0.1, 0.5, 0.5, 5, 5, 0.1 };
                    double LS_LG = 0.7;
                    double Mu_rate = 0.5;
                    Iteration = Convert.ToInt32((O * (1 - LS_LG)) / NPop) + 1;
                    int I_str = Convert.ToInt32(O * LS_LG / Iteration);
                    int best_index = 0;
                    int local_search = 0;
                    double[,] Diversity = new double[Ntest, Iteration];
                    double[,] Fit = new double[Ntest, Iteration];
                    string ss2 = "fit_N_" + No_nonanchor.ToString() + "_Mr_" + Mu_rate.ToString() + "LSLG" + LS_LG.ToString() + "_Pop_" + NPop.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "Struct_" + Struct +range.ToString()+ ".txt";
                    //string ss = "Diver_N_" + No_nonanchor.ToString() + "_Mr_" + Mu_rate.ToString() + "LSLG" + LS_LG.ToString() + "_Pop_" + NPop.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "Struct_" + Struct + "_OF" + OF.ToString() + ".txt";
                    //string ss1 = "Fitness_N_" + No_nonanchor.ToString() + "_Mr_" + Mu_rate.ToString() + "LSLG" + LS_LG.ToString() + "_Pop_" + NPop.ToString() + method + "_Iteration_" + Iteration.ToString() + "_test_" + Ntest.ToString() + "Struct_" + Struct + "_OF" + OF.ToString() + ".txt";
                    if (!File.Exists(ss2))
                    {
                        for (int tt = 0; tt < Ntest; tt++)
                        {
                            indi[] Pop = new indi[NPop];
                            indi Best_nodes = new indi();
                            indi Breed = new indi();
                            Best_nodes.cost = 1E250;
                            for (int ii = 0; ii < NPop; ii++)
                            {
                                Pop[ii] = NFC.Initialization(No_nonanchor, Rand);
                                Pop[ii].cost = NFC.compute_cost(Pop[ii].x, Pop[ii].y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);
                                Pop[ii].local_search = true;

                            }
                            Pop = NFC.mergesort(Pop, NPop);
                            for (int it = 0; it < Iteration; it++)
                            {
                                //perform Steady State GA
                                //for (int ii = 0; ii < NPop; ii++)
                                //{
                                int parent1 = Rand.Next(NPop);
                                int parent2 = -1;
                                if (Struct != "basic")
                                {
                                    int[] S = NFC.Structures(parent1, NPop, Struct);
                                    parent2 = NFC.negative_assortative_mating_strategy2(Pop[parent1].x, Pop[parent1].y, Pop, Rand, S);
                                }
                                else
                                {
                                    parent2 = NFC.negative_assortative_mating_strategy(Pop[parent1].x, Pop[parent1].y, Pop, Rand);
                                }
                                Breed.x = NFC.BLX_Crossover(Pop[parent1].x, Pop[parent2].x);
                                Breed.y = NFC.BLX_Crossover(Pop[parent1].y, Pop[parent2].y);
                                if (Rand.NextDouble() < Mu_rate)
                                    NFC.Mutation_BGA2(ref Breed.x, ref Breed.y, Pop, No_nonanchor, Rand);
                                Breed.cost = NFC.compute_cost(Breed.x, Breed.y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);

                                NFC.Node_treatment_localization(ref Breed.x, ref  Breed.y, No_nonanchor, Rand);
                                Pop = NFC.mergesort(Pop, NPop);
                                if (Breed.cost < Pop[NPop - 1].cost)
                                {
                                    Pop[NPop - 1].x = (double[])Breed.x.Clone();
                                    Pop[NPop - 1].y = (double[])Breed.y.Clone();
                                    Pop[NPop - 1].cost = Breed.cost;
                                    Pop[NPop - 1].local_search = true;
                                }
                                //
                                //Pick the best individual of Population
                                if (local_search == NPop - 1)
                                {
                                    for (int ii = 1; ii < NPop; ii++)
                                    {
                                        Pop[ii] = NFC.Initialization(No_nonanchor, Rand);
                                        Pop[ii].cost = NFC.compute_cost(Pop[ii].x, Pop[ii].y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);
                                        Pop[ii].local_search = true;
                                        Pop[ii].visited = 0;
                                    }
                                    local_search = 0;
                                }

                                indi Best_Indi_BLS = new indi();
                                Pop = NFC.mergesort(Pop, NPop);
                                for (int ii = 0; ii < NPop; ii++)
                                    if (Pop[ii].local_search)
                                    {
                                        Best_Indi_BLS.cost = Pop[ii].cost;
                                        Best_Indi_BLS.x = (double[])Pop[ii].x.Clone();
                                        Best_Indi_BLS.y = (double[])Pop[ii].y.Clone();
                                        Best_Indi_BLS.visited = Pop[ii].visited;
                                        Best_Indi_BLS.rho = Pop[ii].rho;
                                        if (Pop[ii].visited == 1)
                                        {
                                            Best_Indi_BLS.biasx = (double[])Pop[ii].biasx.Clone();
                                            Best_Indi_BLS.biasy = (double[])Pop[ii].biasy.Clone();
                                        }
                                        best_index = ii;
                                        break;
                                    }
                                // Initialize Bias and Rho
                                if (Best_Indi_BLS.visited == 0)
                                {
                                    Best_Indi_BLS.biasx = new double[No_nonanchor];
                                    Best_Indi_BLS.biasy = new double[No_nonanchor];
                                    Best_Indi_BLS.rho = 0.01;
                                }
                                // Local Search Operator
                                indi Best_Indi_ALS = new indi();
                                Best_Indi_ALS.x = (double[])Best_Indi_BLS.x.Clone();
                                Best_Indi_ALS.y = (double[])Best_Indi_BLS.y.Clone();
                                Best_Indi_ALS.cost = Best_Indi_BLS.cost;
                                Best_Indi_ALS.rho = Best_Indi_BLS.rho;
                                Best_Indi_ALS.visited = Best_Indi_BLS.visited;
                                Best_Indi_ALS.biasx = (double[])Best_Indi_BLS.biasx.Clone();
                                Best_Indi_ALS.biasy = (double[])Best_Indi_BLS.biasy.Clone();
                                NFC.Subgrouping_SolisWets_localization(ref Best_Indi_ALS.x, ref Best_Indi_ALS.y, ref Best_Indi_ALS.cost, No_nonanchor, ref Best_Indi_ALS.biasx, ref Best_Indi_ALS.biasy, ref Best_Indi_ALS.rho, I_str, Rand, range, Realx, Realy, gaussian, Noise_Factor, Anchorx, Anchory);
                                // Remove Unchanged particle from population
                                if (Best_Indi_BLS.cost - Best_Indi_ALS.cost == 0)
                                {
                                    Pop[best_index].local_search = false;
                                    //  Pop[best_index].visited = 1;
                                    local_search++;
                                }
                                else
                                {
                                    Pop[best_index].x = (double[])Best_Indi_ALS.x.Clone();
                                    Pop[best_index].y = (double[])Best_Indi_ALS.y.Clone();
                                    Pop[best_index].visited = 1;
                                    Pop[best_index].cost = Best_Indi_ALS.cost;
                                    Pop[best_index].rho = Best_Indi_ALS.rho;
                                    Pop[best_index].biasx = (double[])Best_Indi_ALS.biasx.Clone();
                                    Pop[best_index].biasy = (double[])Best_Indi_ALS.biasy.Clone();
                                    Pop[best_index].local_search = true;
                                }

                                Pop = NFC.mergesort(Pop, NPop);
                                if (Best_nodes.cost > Pop[0].cost)
                                {
                                    Best_nodes.cost = Pop[0].cost;
                                    Best_nodes.x = (double[])Pop[0].x.Clone();
                                    Best_nodes.y = (double[])Pop[0].y.Clone();
                                }

                                /////////////////////
                                //Diversity[tt, it] = NFC.diversity(Pop, NPop, No_nonanchor);
                                //Fit[tt, it] = Best_nodes.cost ;
                            }
                            costs[tt] = NFC.Calc_error(Best_nodes.x, Best_nodes.y, Realx, Realy, range);
                        }


                        sw = new StreamWriter(ss2);
                        for (int tt = 0; tt < Ntest; tt++)
                        {
                            sw.Write(costs[tt].ToString() + " ");
                        }
                        sw.Close();
                        //sw = new StreamWriter(ss);
                        //for (int tt = 0; tt < Ntest; tt++)
                        //{
                        //    for (int ii = 0; ii < Iteration; ii++)
                        //    {
                        //        sw.Write(Diversity[tt, ii].ToString() + " ");
                        //    }
                        //    sw.WriteLine();
                        //}

                        //sw.Close();
                        //sw = new StreamWriter(ss1);
                        //for (int tt = 0; tt < Ntest; tt++)
                        //{
                        //    for (int ii = 0; ii < Iteration; ii++)
                        //    {
                        //        sw.Write(Fit[tt, ii].ToString() + " ");
                        //    }
                        //    sw.WriteLine();
                        //}
                        sw.Close();
                    }
                    //}
                }

            }

            #endregion
            #region FMOA
            if (method == "FMOA")
            {
                double Jr = 0.4;
                double intensity = 1;
                double alpha = 0.6;
                string structure = "Cellular";
                string ss = "Result_Pop_" + NPop.ToString() + ".txt";
                if (!File.Exists(ss))
                {
                    for (int tt = 0; tt < Ntest; tt++)
                    {
                        indi[] Nodes = new indi[NPop];
                        indi[] ONodes = new indi[NPop];
                        indi[] Forces = new indi[NPop];
                        indi[] a = new indi[NPop];
                        indi[] v = new indi[NPop];
                        double[] Magnet = new double[NPop];
                        double[] Mass = new double[NPop];
                        indi Best_indi = new indi();
                        double[] distances = new double[NPop];
                        Best_indi.cost = -1E290;
                        for (int ii = 0; ii < NPop; ii++)
                        {
                            Forces[ii].x = new double[No_nonanchor];
                            a[ii].x = new double[No_nonanchor];
                            v[ii].x = new double[No_nonanchor];
                            Forces[ii].y = new double[No_nonanchor];
                            a[ii].y = new double[No_nonanchor];
                            v[ii].y = new double[No_nonanchor];
                        }
                        for (int ii = 0; ii < NPop; ii++)
                        {
                            Nodes[ii] = NFC.Initialization(No_nonanchor, Rand);
                            ONodes[ii].x = NFC.Opposition_Learning(Nodes[ii].x);
                            ONodes[ii].y = NFC.Opposition_Learning(Nodes[ii].y);
                        }
                        for (int it = 0; it < Iteration; it++)
                        {
                            double Fmax = -1E290; double Fmin = 1000000000;
                            int min_index = 0; int max_index = 0;
                            for (int ii = 0; ii < NPop; ii++)
                            {
                                Nodes[ii].cost = NFC.compute_cost(Nodes[ii].x, Nodes[ii].y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);

                                if (Jr > Rand.NextDouble())
                                {

                                    ONodes[ii].x = NFC.generating_jumping(Nodes, NPop, No_nonanchor, ii);
                                    ONodes[ii].y = NFC.generating_jumping(Nodes, NPop, No_nonanchor, ii);
                                    ONodes[ii].cost = NFC.compute_cost(Nodes[ii].x, Nodes[ii].y, Anchorx, Anchory, range, range, Realx, Realy, Noise_Factor, gaussian);
                                    if (Nodes[ii].cost < ONodes[ii].cost)
                                    {
                                        Nodes[ii].x = (double[])ONodes[ii].x.Clone();
                                        Nodes[ii].y = (double[])ONodes[ii].y.Clone();
                                        Nodes[ii].cost = ONodes[ii].cost;
                                    }
                                    NFC.Node_treatment_localization(ref Nodes[ii].x, ref Nodes[ii].y, No_nonanchor, Rand);
                                    NFC.Node_treatment_localization(ref v[ii].x, ref v[ii].y, No_nonanchor, Rand);
                                }

                                if (Fmax < Nodes[ii].cost)
                                {
                                    Fmax = Nodes[ii].cost;
                                    max_index = ii;
                                }
                                if (Fmin > Nodes[ii].cost)
                                {
                                    Fmin = Nodes[ii].cost;
                                    min_index = ii;
                                }
                            }

                            if (Best_indi.cost < Fmax)
                            {
                                Best_indi.cost = Fmax;
                                Best_indi.x = (double[])Nodes[max_index].x.Clone();
                                Best_indi.y = (double[])Nodes[max_index].x.Clone();
                            }
                            double ranges = Fmax - Fmin;
                            if (ranges == 0)
                            {
                                Fmax = -1000000000; Fmin = 1000000000;
                                min_index = 0; max_index = 0;
                                for (int ii = 0; ii < NPop; ii++)
                                {
                                    Nodes[ii] = NFC.Initialization(No_nonanchor, Rand);
                                    ONodes[ii].x = NFC.Opposition_Learning(Nodes[ii].x);
                                    ONodes[ii].y = NFC.Opposition_Learning(Nodes[ii].y);
                                }
                            }
                            for (int ii = 0; ii < NPop; ii++)
                            {
                                Magnet[ii] = (Nodes[ii].cost - Fmin) / (range);
                                Mass[ii] = Magnet[ii] * intensity + alpha;
                            }
                            for (int ii = 0; ii < NPop; ii++)
                            {
                                Forces[ii].x = new double[No_nonanchor];
                                Forces[ii].y = new double[No_nonanchor];
                            }
                            for (int ii = 0; ii < NPop; ii++)
                            {
                                int[] Neigbour = NFC.Structures(ii, NPop, structure);
                                for (int jj = 0; jj < Neigbour.Length; jj++)
                                {
                                    if (ii != Neigbour[jj])
                                    {

                                        for (int n = 0; n < No_nonanchor; n++)
                                        {
                                            Forces[ii].x[n] = Forces[ii].x[n] + (Nodes[Neigbour[jj]].x[n] - Nodes[ii].x[n]) * Magnet[Neigbour[jj]] / distances[ii];
                                            Forces[ii].y[n] = Forces[ii].y[n] + (Nodes[Neigbour[jj]].y[n] - Nodes[ii].y[n]) * Magnet[Neigbour[jj]] / distances[ii];
                                        }



                                    }
                                }
                            }
                            for (int ii = 0; ii < NPop; ii++)
                            {
                                for (int n = 0; n < No_nonanchor; n++)
                                {
                                    v[ii].x[n] = (Forces[ii].x[n] / Mass[ii]);
                                    v[ii].y[n] = (Forces[ii].y[n] / Mass[ii]);
                                    Nodes[ii].x[n] = Nodes[ii].x[n] + v[ii].x[n];
                                    Nodes[ii].y[n] = Nodes[ii].y[n] + v[ii].y[n];
                                }

                                NFC.Node_treatment_localization(ref Nodes[ii].x, ref Nodes[ii].y, No_nonanchor, Rand);
                                NFC.Node_treatment_localization(ref v[ii].x, ref v[ii].y, No_nonanchor, Rand);
                            }

                        }
                        costs[tt] = Best_indi.cost;
                    }
                    sw = new StreamWriter(ss);
                    for (int tt = 0; tt < Ntest; tt++)
                    {
                        sw.Write(costs[tt].ToString() + " ");
                    }
                    sw.Close();
                }
                //}
                //}
            }
            #endregion
        }


    }
}
