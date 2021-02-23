package comp559.particle;

public class RK4 implements Integrator {
    
    @Override
    public String getName() {
        return "RK4";
    }

    public double[] k1;
    public double[] k2;
    public double[] k3;
    public double[] k4;

    @Override
    public void step(double[] p, int n, double t, double h, double[] pout, Function derivs) {
        // TODO: Objective 6, implement the RK4 integration method

        if ( k1 == null || k1.length != n ) {
            k1 = new double[n]; // beginning of the interval
            k2 = new double[n]; // midpoint of the interval
            k3 = new double[n]; // same but with k2
            k4 = new double[n]; // end of the interval
        }

        derivs.derivs(t, p, k1);

        for(int i = 0; i < n; i++)
        {
            k2[i] = p[i] + h*k1[i]/2.0;
        }

        derivs.derivs(t+h/2.0, k2, k2);

        for(int i = 0; i < n; i++)
        {
            k3[i] = p[i] + h*k2[i]/2.0;
        }

        derivs.derivs(t+(h/2.0), k3, k3);

        for(int i = 0; i < n; i++)
        {
            k4[i] = p[i] + h*k3[i];
        }

        derivs.derivs(t+h, k4, k4);

        for(int i = 0; i < n; i++)
        {
            pout[i] = p[i] + (h/6.0)*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
        }
    }
}
