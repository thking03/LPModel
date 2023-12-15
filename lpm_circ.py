"""
Lumped Parameter Model of the Circulation for BME 307 Final Project
Tyler King, Katy Price

Model architecture, parameters, and initial conditions from Tang et al., 2020: https://link.springer.com/article/10.1007/s12265-020-09953-y. 
    (see GitHub here: https://github.com/tanghongdlut/An-Analog-Circuit-Platform-for-Simulating-Typical-Pulmonary-Hypertensions-from-Point-of-Hemodynamics/blob/master/Entry_for_NormalHumanCirculationSystem.m)
Additional model parameters (coronary and cephalic pathways are fron other papers and integrated into the model).

This script provides code for easy construction of a lumped parameter model of the circulation as described by Tang et al.
Additionally, this model allows for tracking molecule movement through the compartments.

TO-DO:
    - Add in diffusion mode for LPM units
    - Add inductors back into the model
    - Figure out how to extract flow and pressure parameters from the model
"""

import numpy as np
from scipy.integrate import odeint


# Class definitions to easily create models

class LPModel:
    # A dictionary of LPM units of variothius types

    def __init__(self):
        self.units = {}
        pass

    def add_unit(self, unit):
        if not isinstance(unit, LPunit):
            raise Exception("May only add LPunits to the model.")
        self.units[unit.id] = unit

    def get_graph(self):
        # Retrieves a graph of the model for visualization purposes. Only retrieves forward connections.
        if size == 0:
            raise Exception("There are no units in the LPM.")
        graph = []
        unit_dict = self.units # Reference so that we don't need to call multiple times
        for key in unit_dict.keys():
            this_unit = unit_dict[key] # Copy for fast access
            graph.append(this_unit.id, this_unit.after)
        return graph

    def solve_model(self, dt, T):
        # dt is the number of time steps per second.
        # T is the number of seconds to run for
        # We track the volume of each unit in time & flow through inductor if inductor exists.
        unit_dict = self.units

        # Check if model is ready to be run
        break_model = False
        order = sorted(list(unit_dict.keys())) # For now MUST go from 0 to N-1 for N nodes
        for key in order:
            this_unit = unit_dict[key] # Copy for fast access
            if hasattr(this_unit, "name"):
                uid = this_unit.name # Use name for easy reference if exists...
            else:
                uid = this_unit.id # ... or just use ID
            # Check for graph connectedness
            if not hasattr(this_unit, "prev"):
                print("Unit {} has no backward connections.".format(uid))
                break_model = True
            if not hasattr(this_unit, "after"):
                print("Unit {} has no forward connections.".format(uid))
                break_model = True
            # Check initial conditions
            if not hasattr(this_unit, "vstate"):
                print("Unit {} does not have initial volume.".format(uid))
                break_model = True
            # Make sure that all diffusive units have diffusion equations
            if this_unit.Diff and not hasattr(this_unit, "dfunc"):
                print("Unit {} is flagged for diffusion but has no defined equation.".format(uid))
                break_model = True
            # if this_unit.L != 0 and not hasattr(this_unit, "lstate"):
            #     print("Unit {} has inductor with no initial flow.".format(uid))
        if break_model:
            raise Exception("The model is not complete (see printouts).")

        # Run model
        #   - Extract and store dV for each node
        #   - The general form is dV = dt*(SUM[Qin] - SUM[Qout]) = dt*(SUM[(P_prev_i - P_this)/R_prev-i] - SUM[(P_this - P_next_i)/R_this)]
        #   - Obtain P from V using P-V relationships
        #   - For inductor nodes also need to store dQ through inductor w/ form dQ = dt*((P_this - P_next - Q_this*R_this)/L_this)
        #   - We also use the flow terms to do diffusion.
        
        v0 = [] # Initial conditions for volume
        l0 = [] # Initial conditions for inductor flow
        n0 = [] # Initial conditions for diffusion (we track number per compartment and calc. concentration dynamically)
        for key in order:
            v0.append(unit_dict[key].vstate)
            # if unit_dict[key].L != 0:
            #    l0.append(unit_dict[key].lstate)
            n0.append(unit_dict[key].nstate)
        
        def calc_derivative(y, t):
            vols = y[:len(v0)] # Get volumes
            ls = y[len(v0):len(v0)+len(l0)] # Get inductances
            mols = y[len(v0)+len(l0):] # Get current number of molecules
            dvs = []
            dls = []
            dns = []
            vi = 0
            li = 0
            ni = 0
            
            # This does all flow
            for key in order:
                this_unit = unit_dict[key] # access
                this_vol = vols[vi]

                # Get pressure of this unit
                if this_unit.nonlinear[1]:
                    this_P = this_unit.C(this_vol, t=t) # All our PV equations RETURN P
                else:
                    this_P = this_vol/this_unit.C

                # Get the concentration of the unit
                this_conc = mols[ni]/this_vol

                # Get derivatives
                total_flow = 0.0 # total flow is derivative of volume
                total_mass_flow = 0.0 # total mass flow is the derivative of number
                
                if this_unit.Diff:
                    total_mass_flow += this_unit.dfunc(t, this_conc)
                
                for prev_id in this_unit.prev:
                    diode_state = True # Reset the diode state for next calculation
                    # Get pressure of previous unit
                    if unit_dict[prev_id].nonlinear[1]:
                        prev_P = unit_dict[prev_id].C(vols[prev_id], t=t) # All our PV equations RETURN P
                    else:
                        prev_P = vols[prev_id]/unit_dict[prev_id].C
                    # Get the concentration of the previous unit
                    prev_conc = mols[prev_id]/vols[prev_id]
                    # Calculate flow into unit
                    if unit_dict[prev_id].nonlinear[0]:
                        inflow = (prev_P - this_P)/unit_dict[prev_id].R(vols[prev_id], t=t)
                    else:
                        inflow = (prev_P - this_P)/unit_dict[prev_id].R

                    # Solve for the diode state & do bulk fluid flow tracking
                    if unit_dict[prev_id].D_out and t > 0:
                        diode_state = prev_P > this_P # Sets the bool to true if flow goes the right way
                        # print("Got diode state: {} for unit {} at {} mmHg and unit {} at {} mmHg .".format(diode_state, key, this_P, prev_id, prev_P))
                    total_flow += inflow*diode_state # No flow if bool is false bc of valve
                    
                    # Do solute number tracking
                    if inflow >= 0: # This means flow is from previous unit INTO this unit
                        total_mass_flow += inflow*prev_conc*diode_state
                    else:
                        total_mass_flow += inflow*this_conc*diode_state

                    
                for next_id in this_unit.after:
                    diode_state = True # Reset the diode state for next calculation
                    # Get pressure of next unit
                    if unit_dict[next_id].nonlinear[1]:
                        next_P = unit_dict[next_id].C(vols[next_id], t=t) # All our PV equations RETURN P
                    else:
                        next_P = vols[next_id]/unit_dict[next_id].C
                    # Get the concentration of the next unit
                    next_conc = mols[next_id]/vols[next_id]
                    # Calculate flow out of unit
                    if this_unit.nonlinear[0]:
                        outflow = (this_P - next_P)/this_unit.R(this_vol, t=t)
                    else:
                        outflow = (this_P - next_P)/this_unit.R

                    # Solve for the diode state & do bulk fluid flow tracking
                    if this_unit.D_out and t > 0:
                        diode_state = this_P > next_P # Sets the bool to true if flow goes the right way
                        # print("Got diode state: {} for unit {} at {} mmHg and unit {} at {} mmHg .".format(diode_state, key, this_P, next_id, next_P))
                    total_flow -= outflow*diode_state

                    # Do solute number tracking
                    if outflow >= 0: # This means flow is from this unit INTO next unit
                        total_mass_flow -= outflow*this_conc*diode_state
                    else:
                        total_mass_flow -= outflow*next_conc*diode_state
                
                # print("DIAGNOSTIC: for step {}: unit {} has Q = {}".format(t, key, total_flow))
                dvs.append(total_flow)
                dns.append(total_mass_flow)

                vi +=1
                ni +=1

            # print("DIAGNOSTIC: this step has a total volume change of {}".format(sum(dvs)))
            return dvs + dls + dns
        
        t = np.linspace(0, T, int(T/dt)+1)
        states = odeint(calc_derivative, v0 + l0 + n0, t)

        self.last_run = {"states" : states, "times": t} # Save the states for other use
        return states
    
    def get_conc(self):
        if not hasattr(self, "last_run"):
            print("Model has not been run yet. Please run model using self.solve_model().")
        else:
            size = len(self.units.keys())
            vstates = self.last_run["states"][:,:size]
            nstates = self.last_run["states"][:, size:]
            return [nstates[:,i]/vstates[:,i] for i in range(size)]
    
    def get_pressures(self):
        if not hasattr(self, "last_run"):
            print("Model has not been run yet. Please run model using self.solve_model().")
        else:
            size = len(self.units.keys())
            vstates = self.last_run["states"][:,:size]
            times = self.last_run["times"]
            pstates = []
            for key in self.units.keys():
                unit = self.units[key]
                if unit.nonlinear[1]:
                    statevector = []
                    for i, t in enumerate(times):
                        statevector.append(unit.C(vstates[:,key][i], t=t))
                    pstates.append(np.array(statevector)) # This is a vector of pressures
                else:
                    pstates.append(vstates[:,key]/unit.C*np.ones_like(times))
        
        return pstates


class LPunit:
    # A unit in the lumped parameter model

    def __init__(self, id, R, C=0, L=0, D_out=False, links=None, diffusive_region=False, name=None, ics=None):
        # Assign the three standard LINEAR parameters
        self.id = id # use to uniquely identify the unit in the LPM -- maybe could build in some type autonumbering later...?
        self.R = R # There must be vessel resistance
        self.C = C # If this is zero, there is no vessel elastance
        self.L = L # If this is zero, there is no inductance
        self.D_out = D_out # If this is False, there is no valve leaving this unit
        self.Diff = diffusive_region # If this is False, there is NO DIFFUSION into the unit from an outside source
        self.nonlinear = [False, False, False]

        if links:
            self.link_unit(links[0], links[1])
        if name and isinstance(name, str):
            self.name = name
        if ics:
            self.pass_ICs(ics)
    
    def link_unit(self, before, after):
        self.prev = before
        self.after = after

    def pass_ICs(self, ics):
        # Initial conditions should be passed in the form [v_0, n_0] or [v_0, l_0, n_0] if an inductor is present
        # If the system has an inductor it also requires initial inductor flow
        self.vstate = ics[0]
        if self.L != 0:
            self.lstate = ics[1]
            self.nstate = ics[2]
        else:
            self.nstate = ics[1]
    
    def assign_nonlinear(self, nltype, func):
        # Use this to assign NONLINEAR parameters
        if nltype == "R":
            self.R = func
            self.nonlinear[0] = True
        elif nltype == "C":
            self.C = func
            self.nonlinear[1] = True
        elif nltype == "L":
            self.L = func
            self.nonlinear[2] = True

    def assign_diffusion_input(self, func):
        if not self.Diff:
            raise Exception("This unit is not configured to accept a diffusive input. Set self.Diff = True")
        else:
            self.dfunc = func

### When not indicated these formulae and parameters are from Tang et al.
# Key parameters
HR = 60 # Heart rate
[A0, B0, C0] = [0.9, 0.038, 0.145] # Activation functions for the heart
[Aa, Ba, Ca] = [[0.3, 0.35, 0.5, 0.55], [0.045, 0.035, 0.037, 0.036], [0.275, 0.33, 0.375, 0.4]] # Activation functions for the heart
[amin, bmin, ka, kb] = [-2, 0.7, 7, 0.5] # Acivation functions for the heart
[E_esla, E_eslv, E_esra, E_esrv] = [0.3, 4.3, 0.3, 0.8] # Elastances at end-systolic
[M_la, M_lv, M_ra, M_rv] = [0.5, 1.7, 0.5, 0.67] 
[v_0la, v_0lv, v_0ra, v_0rv] = [20, 25, 20, 25]
[v_dla, v_dlv, v_dra, v_drv] = [20, 40, 20, 40]
[l_la, l_lv, l_ra, l_rv] = [0.025, 0.015, 0.025, 0.015]
[Kv, Vsmax] = [40, 3500] # Systemic vein parameters
[n1, n2, k1, k2, kRvc, r0, vc_0, vc_max, vc_min] = [0, -5, 0.15, 0.4, 0.001, 0.025, 130, 350, 50] # Vena cava
[n0, kc, kp1, kp2, krpsa, vsap_min, vsap_max, taop] = [50, 1000, 0.03, 0.2, 0.04, 210, 250, 0.1] # Proximal systemic artery parameters
[f_vaso, f_con] = [0.5, 0.5]

# P-V relationships for the heart using end-diastolic and end-systolic info
# f_con is symphathetic discharge freuqency
# These functions return a pressure given some volume
def etvfunc(t, A, B, C, right=False):
    atv = 0 
    for i in range(4):
        atv += A[i] * np.exp(-(((bmin+kb*f_con)*(t%(60/HR))-(C[i]+right*.013))/B[i])**2)
    return atv

def etafunc(t, A0, B0, C0, right=False):
    return A0 * np.exp(-0.5*((t%(60/HR) - (C0-right*.02))/B0)**2)

def pvt_lv(vol, t):
    pES = (amin+ka*f_con)*E_eslv*(vol-v_dlv)
    pED = M_lv*np.abs(np.exp(l_lv*(vol-v_0lv))-1)
    etvl = etvfunc(t, Aa, Ba, Ca)
    return etvl*pES + (1-etvl)*pED

def pvt_rv(vol, t):
    pES = (amin+ka*f_con)*E_esrv*(vol-v_drv)
    pED = M_rv*np.abs(np.exp(l_rv*(vol-v_0rv))-1)
    etvr = etvfunc(t, Aa, Ba, Ca, right=True)
    return etvr*pES + (1-etvr)*pED

def pvt_la(vol, t):
    pES = E_esla*(vol-v_dla)
    pED = M_la*np.abs(np.exp(l_la*(vol-v_0la))-1)
    etal = etafunc(t, A0, B0, C0)
    return etal*pES + (1-etal)*pED

def pvt_ra(vol, t):
    pES = E_esra*(vol-v_dra)
    pED = M_ra*np.abs(np.exp(l_ra*(vol-v_0ra))-1)
    etar = etafunc(t, A0, B0, C0, right=True)
    return etar*pES + (1-etar)*pED

# P-V relationships for some vessels
# Write everything as a function of volume and time so time is accepted
def pv_sysv(vol, t):
    return -Kv*np.log10(Vsmax/vol - 0.99)

def pv_vena(vol, t):
    if vol >= vc_0:
        return n1 + k1*(vol-vc_0)
    else:
        return n2 + k2*np.exp(vol/vc_min)

def rv_vena(vol, t):
    return kRvc*(vc_max/vol)+r0

# Proximal systemic artery is complex...
# f_vaso is vasoconstriction (sympathetic frequency)
def pv_proxart(vol, t):
    pv_proxart_a = kc*np.log10((vol-vsap_min)/n0+1)
    pv_proxart_p = kp1*np.exp(taop*(vol-vsap_min)) + kp2*(vol-vsap_min)**2
    return f_vaso*pv_proxart_a + (1-f_vaso)*pv_proxart_p
# ... also has resistance relationship
def rv_proxart(vol, t):
    return krpsa*(np.exp(4*f_vaso)+(vsap_max/vol)**2)