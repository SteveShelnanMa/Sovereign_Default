
################################################################
##                   University of Rochester                  ## 
##                        Xiaonan Ma                          ##
##                       19 Jan, 2021                         ##
################################################################



!pip install --upgrade quantecon
import matplotlib.pyplot as plt
import mpl_toolkits.axisartist as axisartist
import numpy as np
import quantecon as qe
import random
from numba import jit, jitclass, int64, float64
%matplotlib inline

#################################################################
#  In this part:                                                #
#     setting the parameter values up,                          #
#     setting up the value functions in terms of default,       #
#     discretizing the AR(1) process for y using Tauchen method #
#     and setting up the grids for state variables B and y      #
#################################################################

# Define the data information for the jitclass
arellano_data = [
    ('B', float64[:]), ('P', float64[:, :]), ('y', float64[:]),
    ('β', float64), ('σ', float64), ('r', float64),
    ('ρ', float64), ('η', float64), ('θ', float64),
    ('def_y', float64[:]),
    ('y_hat',float64)
]

# Define utility function
@jit(nopython=True)
def u(c, σ):
    return c**(1-σ)/(1-σ)

@jitclass(arellano_data)
class Arellano_Economy:
    """
    Parameters:
    B : vector(float64)
        A grid for bond holdings
    P : matrix(float64)
        The transition matrix for a country's output
    y : vector(float64)
        The possible output states
    β : float
        Discount factor, targeting 3% default probability
    σ : float
        Risk-aversion
    r : float
        Risk free interest rate(U.S. 5 year bond quarterly yield)
    θ : float
        Probability of re-entering financial markets in each period
    y_hat: float
        Output costs, targeting 5.53% debt service to GDP

    Stochastic structure based on Argentina’s GDP
    ρ : float
        Persistence in the income process
    η : float
        Standard deviation of the income process
    """

    def __init__(
            self, B, P, y,
            β=0.953, σ=2.0, r=0.017,
            ρ=0.945, η=0.025, θ=0.282,
            y_hat=0.969, figure2=False
        ):
        # Save parameters   
            self.B, self.P, self.y = B, P, y
            self.β, self.σ, self.r = β, σ, r
            self.ρ, self.η, self.θ = ρ, η, θ
            self.y_hat = y_hat
            self.def_y = np.minimum(self.y_hat * np.mean(y), y)

    def bellman_default(self, iy, EVd, EV):
        """
        Return the value of the Bellman equation when the country is 
        in a defaulted state on their debt,
        i.e. equ.(8) in Arellano(2008)
        """
        # Unpack certain parameters for simplification
        β, σ, θ = self.β, self.σ, self.θ

        # Compute continuation value
        zero_ind = len(self.B) // 2
        cont_value = θ * EV[iy, zero_ind] + (1 - θ) * EVd[iy]

        return u(self.def_y[iy], σ) + β*cont_value

    def bellman_nondefault(self, iy, iB, q, EV, iB_tp1_star=-1):
        """
        Return the value of the Bellman equation when the country is 
        not in a defaulted state,
        i.e. equ.(9) in Arellano(2008)
        """
        # Unpack certain parameters for simplification
        β, σ, θ = self.β, self.σ, self.θ
        B, y = self.B, self.y

        # Compute the RHS of Bellman equation
        if iB_tp1_star < 0:
            iB_tp1_star = self.compute_savings_policy(iy, iB, q, EV)
        c = max(y[iy] - q[iy, iB_tp1_star]*B[iB_tp1_star] + B[iB], 1e-14)

        return u(c, σ) + β*EV[iy, iB_tp1_star]

    def compute_savings_policy(self, iy, iB, q, EV):
        """
        Finds the index of the debt/savings that maximizes the value 
        function for a particular state given prices and value function
        """
        # Unpack certain parameters for simplification
        β, σ, θ = self.β, self.σ, self.θ
        B, y = self.B, self.y

        # Compute the RHS of Bellman equation
        current_max = -1e14
        iB_tp1_star = 0
        for iB_tp1, B_tp1 in enumerate(B):
            c = max(y[iy] - q[iy, iB_tp1]*B[iB_tp1] + B[iB], 1e-14)
            m = u(c, σ) + β*EV[iy, iB_tp1]

            if m > current_max:
                iB_tp1_star = iB_tp1
                current_max = m

        return iB_tp1_star



β, σ, r = 0.953, 2.0, 0.017
ρ, η, θ1, θ2 = 0.945, 0.025, 0.282, 0
ny = 21
nB = 200

y_hat1, y_hat2 = 0.969, float('inf')
B_lowbar, B_upbar = -2.5, 0.2  

# seting up the grids for state variables B
Bgrid = np.linspace(B_lowbar, B_upbar, nB)

# discretizing the AR(1) process for y
mc = qe.markov.tauchen(ρ, η, 0, 3, ny)
# seting up the grids for state variables y (log-normal)
# P is the transition matrix drived from tauchen method
ygrid = np.exp(mc.state_values)
P = mc.P
P2 = np.ones((ny,ny))/ny

# parameters for figure 3 & 4
ae = Arellano_Economy(
    Bgrid, P, ygrid, β=β, σ=σ, r=r, ρ=ρ, η=η, θ=θ1,
    y_hat = y_hat1
)

# parameters for figure 2: y is iid, θ=0, h(y)=y
ae2 = Arellano_Economy(
    Bgrid, P2, ygrid, β=β, σ=σ, r=r, ρ=ρ, η=η, θ=θ2,
    y_hat = y_hat2
)

#################################################################
#  In this part:                                                #
#     setting up the initial guesses in the solve() function,   #
#     solving the model using VFI,                              #
#     and replicating figures                                   #
#################################################################

# solve the model using VFI
@jit(nopython=True)
def solve(model, tol=1e-05, maxiter=10_000):
    """
    Given an Arellano_Economy type, this function computes the optimal
    policy and value functions
    """
    # Unpack certain parameters for simplification
    β, σ, r, θ = model.β, model.σ, model.r, model.θ
    # np.ascontinguousarray() returns a contiguous array in memory(C order)
    # to make it run faster
    B = np.ascontiguousarray(model.B)
    P, y = np.ascontiguousarray(model.P), np.ascontiguousarray(model.y)
    nB, ny = B.size, y.size

    # Allocate space and set up the initial guesses
    iBstar = np.zeros((ny, nB), int64)
    default_prob = np.zeros((ny, nB))
    default_states = np.zeros((ny, nB))
    q = np.ones((ny, nB)) * 0.95
    Vd = np.zeros(ny)
    Vc, V, Vupd = np.zeros((ny, nB)), np.zeros((ny, nB)), np.zeros((ny, nB))

    it = 0
    dist = 10.0
    while (it < maxiter) and (dist > tol):

        # Compute expectations used for this iteration
        EV = P@V
        EVd = P@Vd

        for iy in range(ny):
            # Update value function for default state
            Vd[iy] = model.bellman_default(iy, EVd, EV)

            for iB in range(nB):
                # Update value function for non-default state
                iBstar[iy, iB] = model.compute_savings_policy(iy, iB, q, EV)
                Vc[iy, iB] = model.bellman_nondefault(iy, iB, q, EV, iBstar[iy, iB])

        # Once value functions are updated, can combine them to get
        # the full value function
        Vd_compat = np.reshape(np.repeat(Vd, nB), (ny, nB))
        Vupd[:, :] = np.maximum(Vc, Vd_compat)

        # Can also compute default states and update prices
        default_states[:, :] = 1.0 * (Vd_compat > Vc)
        default_prob[:, :] = P @ default_states
        q[:, :] = (1 - default_prob) / (1 + r)

        # Check tolerance etc...
        dist = np.max((np.abs(Vupd - V))/(1+np.abs(V)))
        V[:, :] = Vupd[:, :]
        it += 1
        if it % 50==0:
            print("Current iteration number: ",it)
    print("Total iteration number: ", it)

    return V, Vc, Vd, iBstar, default_prob, default_states, q

# return the results and running time
import datetime
start = datetime.datetime.now()

V, Vc, Vd, iBstar, default_prob, default_states, q = solve(ae)
V2, Vc2, Vd2, iBstar2, default_prob2, default_states2, q2 = solve(ae2)

end = datetime.datetime.now()
print ("Running time:" ,(end-start))


# now replicate figures 

# [figure 2]
qB = np.transpose(np.reshape(np.repeat(Bgrid, ny),(nB,ny))) * q2

# get the q(B)*B and denote it as q_B
def find_q_and_lowBi(qB):
    q_B = np.zeros(qB.shape[1])   
    for i in range(qB.shape[1]):
        q_B[i]=qB[:,i].mean()
    return q_B

q_B = find_q_and_lowBi(qB)

# find the index of B for the lowest element of q(B)*B 
qBmin, lowB_i = min((val, idx) for (idx, val) in enumerate(q_B))

fig = plt.figure(figsize=(9, 6))
ax = axisartist.Subplot(fig, 111)  
fig.add_axes(ax)

ax.set_title("Figure 2: Total Resources Borrowed")
ax=plt.gca()  # gca:get current axis
# Sets the right and top borders of the image to not be displayed
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
# move the locations of x&y-axis
ax.spines['bottom'].set_position(('data',0))  
ax.spines['left'].set_position(('data',0)) 

ax.axis[:].set_visible(False)
ax.axis["x"] = ax.new_floating_axis(0,0)
ax.axis["x"].set_axisline_style("->", size = 1.0)
ax.axis["y"] = ax.new_floating_axis(1,0)
ax.axis["y"].set_axisline_style("->", size = 1.0)

ax.axis["x"].set_axis_direction("top")
ax.axis["y"].set_axis_direction("right")
ax.axis["y"].set_axis_direction("left")

ax.annotate("", xy=(0.4, 0.4), xytext=(-0.75, -0.4),# unnecessary, just to make the plot prettier
    arrowprops=dict(arrowstyle="->",connectionstyle="angle3,angleA=0,angleB=-90"))
ax.annotate(r"$slope=\frac{1}{1+r}$", xy=(-0.41, -0.4), xytext=(-1.25, -0.4),
        arrowprops=dict(arrowstyle="->",connectionstyle="angle3,angleA=0,angleB=-90"))
ax.annotate('B*',xy=(Bgrid[lowB_i]-0.015, 0.01), xytext=(-.62, .06))
ax.annotate("$q(B')B'$",xy=(0,.2),xytext=(0.05,.3))
ax.annotate("$B'$",xy=(0,.2),xytext=(-2.25,-.05))

plt.vlines(Bgrid[lowB_i],-0.6,0.3,linestyles = "dashed")
plt.plot([Bgrid[lowB_i],0],[qBmin,qBmin],'--',lw=2)
plt.plot(Bgrid,q_B)
plt.show()




# [figure 3]
# Create "Y High" and "Y Low" values as 5% devs from mean
high, low = np.mean(ae.y) * 1.05, np.mean(ae.y) * .95
iy_high, iy_low = (np.searchsorted(ae.y, x) for x in (high, low))

# Extract a suitable plot grid
x = []
q_low = []
q_high = []

for i in range(nB):
    b = ae.B[i] 
    x.append(b)
    q_low.append(q[iy_low, i])
    q_high.append(q[iy_high, i])

qhigh = np.asarray(q_high)
qlow  = np.asarray(q_low)

'''    display two figures separately:
fig1 = plt.figure(figsize=(9,6))
fig2 = plt.figure(figsize=(9,6))
ax0 = fig1.add_subplot(111)
ax1 = fig2.add_subplot(111)
#fig, (ax0,ax1) = plt.subplots(1,2, figsize=(15,10))
#fig.suptitle("Figure 3: Bond Prices and Assets ")
ax0.plot(x, q_high, label="$y_H$", lw=2, alpha=0.7)
ax0.plot(x, q_low, label="$y_L$", lw=2, alpha=0.7)
ax0.set_xlabel("$B'$")
ax0.legend(loc='upper left', frameon=False)

ax1.plot(ae.B, 1/qlow[iBstar[iy_low]]-1, label="$y_H$", lw=2, alpha=0.7)
ax1.plot(ae.B, 1/qhigh[iBstar[iy_high]]-1, label="$y_L$", lw=2, alpha=0.7)
ax1.set_xlim(-1.5,.2)
ax1.set_ylim(0,.2)
ax1.set_xlabel("$B$")
ax1.legend(loc='upper right', frameon=False)

plt.show(1,2)
'''

# display two figures in a row
width=18
height=6
rows = 1
cols = 2
fig=plt.figure(figsize=(width,height))
axes = []
for i in range(cols*rows): 
    axes.append(fig.add_subplot(rows, cols, i+1) ) 

axes[0].plot(x, q_high, label="$y_H$", lw=2, alpha=0.7)
axes[0].plot(x, q_low, label="$y_L$", lw=2, alpha=0.7)
axes[0].set_xlabel("$B'$")
axes[0].legend(loc='upper left', frameon=False)
axes[0].set_title("Bond Price Schedule q(B',y)")

axes[1].plot(ae.B, 1/qlow[iBstar[iy_low]]-1, label="$y_H$", lw=2, alpha=0.7)
axes[1].plot(ae.B, 1/qhigh[iBstar[iy_high]]-1, label="$y_L$", lw=2, alpha=0.7)
axes[1].set_xlim(-1.5,.25)
axes[1].set_ylim(0.01,.125)
axes[1].set_xlabel("$B$")
axes[1].legend(loc='upper right', frameon=False)
axes[1].set_title("Equilibrium Interest Rate "r"$\frac{1}{q(B'(B),y)}-1$")

#fig.tight_layout()
fig.suptitle("Figure 3: Bond Prices and Assets")
plt.show()




# [figure 4]
# Create "Y High" and "Y Low" values as 5% devs from mean
high, low = np.mean(ae.y) * 1.05, np.mean(ae.y) * .95
iy_high, iy_low = (np.searchsorted(ae.y, x) for x in (high, low))


'''    display two figures separately:
fig, ax = plt.subplots(figsize=(10, 6.5))
ax.plot(ae.B, Bgrid[iBstar[iy_high]], label="$y_H$", lw=2, alpha=0.7)
ax.plot(ae.B, Bgrid[iBstar[iy_low]], label="$y_L$", lw=2, alpha=0.7)
ax.legend(loc='upper left')
ax.set(xlabel="$B$", ylabel=r"$B'(y, B)$")
ax.set_xlim(ae.B.min(), ae.B.max())
plt.show()

fig, ax = plt.subplots(figsize=(10, 6.5))
ax.set_title("Value Functions")
ax.plot(ae.B, V[iy_high], label="$y_H$", lw=2, alpha=0.7)
ax.plot(ae.B, V[iy_low], label="$y_L$", lw=2, alpha=0.7)
ax.legend(loc='upper left')
ax.set(xlabel="$B$", ylabel="$V(y, B)$")
ax.set_xlim(ae.B.min(), ae.B.max())
plt.show()
'''


# display two figures in a row
fig2=plt.figure(figsize=(width,height))
axes2 = []
for i in range(cols*rows): 
    axes2.append(fig2.add_subplot(rows, cols, i+1) ) 

axes2[0].plot(ae.B, Bgrid[iBstar[iy_high]], label="$y_H$", lw=2, alpha=0.7)
axes2[0].plot(ae.B, Bgrid[iBstar[iy_low]], label="$y_L$", lw=2, alpha=0.7)
axes2[0].set(xlabel="$B$", ylabel=r"$B'(y, B)$")
axes2[0].legend(loc='upper left', frameon=False)
axes2[0].set_xlim(ae.B.min(), ae.B.max())
axes2[0].set_title("Savings Function B'(B,y)")

axes2[1].plot(ae.B, V[iy_high], label="$y_H$", lw=2, alpha=0.7)
axes2[1].plot(ae.B, V[iy_low], label="$y_L$", lw=2, alpha=0.7)
axes2[1].legend(loc='upper left', frameon=False)
axes2[1].set_xlim(ae.B.min(), ae.B.max())
axes2[1].set(xlabel="$B$", ylabel=r"$V(B,y)$")
axes2[1].set_title("Value Function V(B,y)")

#fig.tight_layout()
fig2.suptitle("Figure 4: Savings and Value Functions")
plt.show()
