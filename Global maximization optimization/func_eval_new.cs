using System;
using MathNet.Numerics;
using System.Collections;
public class function_evaluation2013
{
    public function_evaluation2013()
	{
	}
    public double f1(int x,int D,double[]data_f0,int initial_flag)
    {

    idx = checkBounds(x, lb, ub);
    x = x-repmat(xopt, 1, ps);
    fit = elliptic(x);
    fit(idx) = NaN;
    }
    public double elliptic(double[]x,int D,double condition)
    {
       
        coefficients = condition .^ linspace(0, 1, D); 
        fit = coefficients * T_irreg(x).^2; 
    }
    public double[] T_irreg(double []f)
    {
        a = 0.1;
        double[]g =(double[])f.Clone(); 
        ArrayList idx=new ArrayList();
        for(int uu=0;uu<f.Length;uu++)
        {
            if(f[uu]>0)
                idx.Add(uu);
        }
        for(int ii=0;ii<idx.Count;ii++)
        {
            g[idx[ii]]= Math.Log(f[idx[ii]])/a;
            g[idx[ii]] =Math.Pow(Math.Exp(g[idx[ii]] + 0.49*(Math.Sin(g[idx[ii]]) +Math.Sin(0.79*g[idx[ii]]))),a);
        }
        idx.Clear();
        for(int uu=0;uu<f.Length;uu++)
        {
            if(f[uu]<0)
                idx.Add(uu);
        }
        for(int ii=0;ii<idx.Count;ii++)
        {
            g[idx[ii]]= Math.Log(-f[idx[ii]])/a;
            g[idx[ii]] =-Math.Pow(Math.Exp(g[idx[ii]] + 0.49*(Math.Sin(g[idx[ii]]) +Math.Sin(0.79*g[idx[ii]]))),a);
        }
        return g;
    }
}
