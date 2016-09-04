import sympy as sym
import numpy as np
import sys, getopt
import argparse
import math as ma

# SCRIPT USAGE:
# This script takes in command line arguments
# V,I,T,w,dt,a,b,u which if not given defaults to certain values
# u is either linear, quadratic or cubic
# and the source term is manufactured from this.
# a and b are additional constants for quadratic and cubic case

V, t, I, w, dt = sym.symbols('V t I w dt')  # global symbols
a, b  = sym.symbols('a b')
f = None # global variable for the source term in the ODE


def ode_source_term(u):
    """Return the terms in the ODE that the source term
    must balance, here u'' + w**2*u.
    u is symbolic Python function of t."""
    return sym.diff(u(t), t, t) + w**2*u(t)

def residual_discrete_eq(u):
    """Return the residual of the discrete eq. with u inserted."""
    R = ode_source_term(u)-DtDt(u,dt)-w**2*u(t)
    return sym.simplify(R)

def residual_discrete_eq_step1(u):
    # R = residual_discrete_eq(u).subs(t,0) <- assume u(-dt) defined
    # using centered difference for u'(0) to obtain u(-dt)
    # gives u(-dt) = u(dt)-2*V*dt, with second order error
    R = ode_source_term(u).subs(t,0)-w**2*u(0)-(2*u(dt)-2*V*dt-2*I)/(dt**2)
    return sym.simplify(R)

def DtDt(u, dt):
    """Return 2nd-order finite difference for u_tt.
    u is a symbolic Python function of t.
    """
    return (u(t+dt)+u(t-dt)-2*u(t))*1/(dt**2)


def main(u):
    """
    Given some chosen solution u (as a function of t, implemented
    as a Python function), use the method of manufactured solutions
    to compute the source term f, and check if u also solves
    the discrete equations.
    """
    print '=== Testing exact solution: %s ===' %u_type
    print "Initial conditions u(0)=%s, u'(0)=%s:" %\
        (u(t).subs(t, 0), sym.diff(u(t), t).subs(t, 0))

    # Method of manufactured solution requires fitting f
    global f  # source term in the ODE
    f = sym.simplify(ode_source_term(u))

    # Residual in discrete equations (should be 0)
    print 'residual step1:', residual_discrete_eq_step1(u)
    print 'residual:', residual_discrete_eq(u)
    print 'Numerical Results:'
    [u_h,t_h,u_err] = solver(V_h,I_h,T_h,w_h,dt_h,a_h,b_h,f,u,u_type)
    print("%s") %np.sqrt(dt_h*np.dot(u_err,u_err))
    #nosetest test_quadratic()
    #

def linear():
    main(lambda t: V*t + I)

def quadratic():
    main(lambda t: a*t**2+V*t+I)

def cubic():
    main(lambda t: a*t**3+b*t**2+V*t+I)

def solver(V_h,I_h,T_h,w_h,dt_h,a_h,b_h,f,u,u_type):
    dt_h = float(dt_h)
    Nt = int(round(T_h/dt_h))
    u_h = np.zeros(Nt+1)
    t_h = np.linspace(0,Nt*dt_h,Nt+1)
    f_lam = sym.lambdify(t,f.subs([(V,V_h),(I,I_h),(w,w_h),(dt,dt_h),(a,a_h),(b,b_h)]),"numpy")
    f_h = f_lam(t_h)
    if u_type == "linear":
        u_e = V_h*t_h+I_h
    elif u_type == "cubic":
        u_e = a_h*np.power(t_h,3)+b_h*np.square(t_h)+V_h*t_h+I_h
    elif u_type == "quadratic":
        u_e = a_h*np.square(t_h)+V_h*t_h+I_h

    u_h[0] = I_h
    u_h[1] = (f_h[0]*dt_h**2 + (2-(w_h*dt_h)**2) * I_h + 2*V_h*dt_h)/2
    for n in range(1,Nt):
        u_h[n+1]= (f_h[n]-w_h**2*u_h[n])*dt_h**2+2*u_h[n]-u_h[n-1]
    u_err = u_h - u_e
    return u_h,t_h,u_err

def test_quadratic():
    u = lambda t: a*t**2+V*t+I
    f = sym.simplify(ode_source_term(u))
    [V_h,I_h,T_h,w_h,dt_h,u_type,a_h,b_h] = [1.0,1.0,1.0,1.0,0.001,"quadratic",1.0,1.0]
    [u_h,t_h,u_err] = solver(V_h,I_h,T_h,w_h,dt_h,a_h,b_h,f,u,u_type)
    assert np.sqrt(dt_h*np.dot(u_err,u_err))<10**-10



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-V", type=float, default=1.0)
    parser.add_argument("-I",type=float, default=1.0)
    parser.add_argument("-T",type=float, default=1.0)
    parser.add_argument("-w",type=float, default=1.0)
    parser.add_argument("-dt",type=float, default=0.5)
    parser.add_argument("-u",type=str, default="linear")
    parser.add_argument("-a",type=float,default=1.0)
    parser.add_argument("-b",type=float,default=1.0)
    args = parser.parse_args()
    V_h = args.V
    I_h = args.I
    T_h = args.T
    w_h = args.w
    dt_h = args.dt
    u_type = args.u
    a_h = args.a
    b_h = args.b
    if u_type=="quadratic":
        quadratic()
    elif u_type=="cubic":
        cubic()
    else:
        linear()
    
    



    