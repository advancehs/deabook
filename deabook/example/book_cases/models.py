"""Small independent LP implementations used for transparent book calculations.

HiGHS solutions are checked for primal feasibility. No observation is silently
removed when a cross-period program is infeasible: that distance becomes NaN.
"""
import numpy as np
from scipy.optimize import linprog


def solve(c, A=None, b=None, E=None, e=None, bounds=(0, None), allow_infeasible=False):
    r = linprog(c, A_ub=A, b_ub=b, A_eq=E, b_eq=e, bounds=bounds, method='highs')
    if not r.success:
        if allow_infeasible and r.status == 2:
            return None
        raise RuntimeError(r.message)
    if A is not None and np.max(np.asarray(A) @ r.x - b) > 2e-6:
        raise ArithmeticError('Inequality residual exceeds tolerance')
    if E is not None and np.max(np.abs(np.asarray(E) @ r.x - e)) > 2e-6:
        raise ArithmeticError('Equality residual exceeds tolerance')
    return r


def radial(X, Y, x, y, rts='crs', orient='input', energy=None, allow_infeasible=False):
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    n, m = X.shape
    if orient == 'input':
        direction = x.copy() if energy is None else np.eye(m)[energy] * x
        A = np.vstack([np.column_stack([X.T, -direction]), np.column_stack([-Y.T, np.zeros(Y.shape[1])])])
        b = np.r_[x - direction, -y]
        c = np.r_[np.zeros(n), 1]
    else:
        A = np.vstack([np.column_stack([X.T, np.zeros(m)]), np.column_stack([-Y.T, y])])
        b = np.r_[x, np.zeros(Y.shape[1])]
        c = np.r_[np.zeros(n), -1]
    E = np.array([np.r_[np.ones(n), 0]]) if rts == 'vrs' else None
    r = solve(c, A, b, E, [1] if E is not None else None, allow_infeasible=allow_infeasible)
    if r is None: return np.nan
    return r.x[-1] if orient == 'input' else 1 / r.x[-1]


def bad_distance(X, Y, B, x, y, b, rts='crs', strong=False, allow_infeasible=False):
    """Bad-output-only proportional DDF; CRS or common-abatement VRS.

    VRS uses the book's common-abatement factor: X'lambda <= theta*x,
    sum(lambda)=theta, 0<=theta<=1. Eliminating theta gives the rows below.
    """
    n = len(X)
    A = np.vstack([np.column_stack([X.T, np.zeros(X.shape[1])]), np.column_stack([-Y.T, np.zeros(Y.shape[1])])])
    rhs = np.r_[x, -y]
    E = np.column_stack([B.T, b])
    if strong:
        A = np.vstack([A, E]); rhs = np.r_[rhs, b]; E = None; eq = None
    else: eq = b
    if rts == 'vrs':
        A[:X.shape[1], :n] -= x[:, None]
        rhs[:X.shape[1]] = 0
        A = np.vstack([A, np.r_[np.ones(n), 0]]); rhs = np.r_[rhs, 1]
    r = solve(np.r_[np.zeros(n), -1], A, rhs, E, eq,
              bounds=[(0, None)]*n+[(None, 1)], allow_infeasible=allow_infeasible)
    return (np.nan, None) if r is None else (1-r.x[-1], r)


def sbm(X, Y, x, y, rts='vrs', B=None, b=None):
    """Tone non-oriented SBM; bad outputs use B lambda + slack_b = t*b."""
    n,m=X.shape; s=Y.shape[1]; k=0 if B is None else B.shape[1]
    c=np.r_[np.zeros(n),-1/(m*x),np.zeros(s+k),1]
    E=np.vstack([
        np.column_stack([X.T,np.eye(m),np.zeros((m,s+k)),-x]),
        np.column_stack([Y.T,np.zeros((s,m)),-np.eye(s),np.zeros((s,k)),-y])])
    if k:
        E=np.vstack([E,np.column_stack([B.T,np.zeros((k,m+s)),np.eye(k),-b])])
    norm=np.r_[np.zeros(n+m),1/((s+k)*y),np.array([]) if k==0 else 1/((s+k)*b),1]
    E=np.vstack([E,norm]);e=np.r_[np.zeros(m+s+k),1]
    if rts=='vrs':E=np.vstack([E,np.r_[np.ones(n),np.zeros(m+s+k),-1]]);e=np.r_[e,0]
    r=solve(c,E=E,e=e)
    t=r.x[-1]
    return r.fun,r.x[n:n+m+s+k]/t,r


def material_balance(X,Y,B,x,y,b,coef,level=1,keep_output=False):
    """Book chapter 9 specialization: nonmaterial K,L, material F, electricity.

    Eliminates nonnegative input/output slacks from the printed MB program.
    Variables: lambda, target emission, optional target F/K/L/Y.
    """
    n=len(X);N=n+1+(level>=2)+(2 if level>=3 else 0)+(level>=4)
    obj=n;fp=n+1 if level>=2 else None;kn=n+2 if level>=3 else None;yp=N-1 if level>=4 else None
    rows=[];rhs=[]
    for j in range(3):
        row=np.zeros(N);row[:n]=X[:,j]
        ix=fp if j==2 else None if kn is None else kn+j
        if ix is None:rhs.append(x[j])
        else:row[ix]=-1;rhs.append(0)
        rows.append(row)
    row=np.zeros(N);row[:n]=-Y[:,0]
    if yp is None:rhs.append(-y[0])
    else:row[yp]=1;rhs.append(0)
    rows.append(row)
    if keep_output and yp is not None:
        row=np.zeros(N);row[yp]=-1;rows.append(row);rhs.append(-y[0])
    row=np.zeros(N);row[:n]=B[:,0];row[obj]=-1;rows.append(row);rhs.append(0)
    eq=np.zeros(N);eq[obj]=1
    if fp is None:eq[:n]=-coef*X[:,2]
    else:eq[fp]=-coef
    r=solve(np.eye(N)[obj],np.array(rows),np.array(rhs),np.array([eq]),np.array([b[0]-coef*x[2]]))
    return r.fun
