#### Funciones de Gating (McCormick & Huguenard, 1992)
#### Gating Functions (McCormick & Huguenard, 1992)

# Activación de Corriente de Sodio / Sodium Current Activation
alpha_m(V::Float64) = 0.091*(V+38)/(1-exp(-(V+38)/5))
beta_m(V::Float64) = -0.062*(V+38)/(1-exp((V+38)/5))
m_inf(V::Float64) = alpha_m(V) / (alpha_m(V) + beta_m(V))
tau_m(V::Float64) = 1 / (alpha_m(V) + beta_m(V))

# Inactivación de Corriente Sodio / Sodium Current Inactivation
alpha_h(V::Float64) = 0.016*exp((-55-V)/15)
beta_h(V::Float64) = 2.07/(exp((17-V)/21)+1)
h_inf(V::Float64) = alpha_h(V) / (alpha_h(V) + beta_h(V))
tau_h(V::Float64) = 1 / (alpha_h(V) + beta_h(V))

# Corriente NaP / Persistent Sodium Current (NaP)
alpha_mp(V::Float64) = 0.091*(V+38)/(1-exp(-(V+38)/5))
beta_mp(V::Float64) = -0.062*(V+38)/(1-exp((V+38)/5))
mp_inf(V::Float64) = 1/(1+exp((-49-V)/5)) 
tau_mp(V::Float64) = 1 / (alpha_mp(V) + beta_mp(V))

# Activación de la corriente A activation / Activation of A-type Current
mA_inf(V::Float64) = 1/(1+exp(-(V+60)/8.5))
tau_mA(V::Float64) = 0.37 + 1/(exp((V+35.82)/19.697)+exp((V+79.69)/-12.7))

# Inactivación de la corriente A / Inactivation of A-type Current
hA_inf(V::Float64) = 1/(1+exp((V+78)/6))
function tau_hA(V::Float64)
    if V < -63
        tau_hA = 1/(exp((V+46.05)/5)+exp((V+238.4)/-37.45))
    else
        tau_hA = 19
    end
    return tau_hA
end

# Activación de la corriente de Calcio tipo T / Activation of T-type Calcium Current
mt_inf(V::Float64) = 1/(1+exp(-(V+57)/6.2))
tau_mt(V::Float64) = 0.612 + 1/(exp(-(V+131.6)/16.7)+exp((V+16.8)/18.2))

# Inactivación de la corriente de Calcio tipo T / Inactivation of T-type Calcium Current
ht_inf(V::Float64) = 1/(1+exp((V+81)/4.03))
function tau_ht(V::Float64)
    if V < -80
        tau_ht = exp((V+467)/66.6)
    else
        tau_ht = exp(-(V+21.88)/10.2)+28
    end
    return tau_ht
end

# Activación de la corriente de Potasio / Activation of Potassium Current
mK2_inf(V::Float64) = 1/(1+exp((V+43)/-17))
tau_mK2(V::Float64) = 1/(exp((V-81)/25.6)+exp((V+132)/-18)) + 9.9

#Inactivación de Potasio tipo A / Inactivation of A-type Potassium Current
hK2a_inf(V::Float64) = 1/(1+exp((V+58)/10.6))
tau_hK2a(V::Float64) = 1/(exp((V-1.329)/200)+exp((V+130)/-7.1))+120

#Inactivación de Potasio tipo B / Inactivation of B-type Potassium Current
hK2b_inf(V::Float64) = 1/(1+exp((V+58)/10.6))
function tau_hK2b(V::Float64)
    if V < -70
        tau_hK2b = 1/(exp((V-1.329)/200)+exp((V+130)/-7.1))+120
    else
        tau_hK2b = 8.9
    end
    return tau_hK2b
end

# Corriente de Potasio tipo-C (activada por Calcio) / C-type Potassium Current (Calcium-activated)
alpha_mc(V::Float64, CaL::Float64) = 2500*CaL*exp(V/24)#2.5e5*CaL*exp(V/24)
beta_mc(V::Float64) = 0.1*exp(-V/24)
mc_inf(V::Float64, CaL::Float64) = alpha_mc(V, CaL) / (alpha_mc(V, CaL) + beta_mc(V))
tau_mc(V::Float64, CaL::Float64) = 1 / (alpha_mc(V, CaL) + beta_mc(V))

# Corriente de Calcio de tipo-L / L-type Calcium Current
alpha_ml(V::Float64) = 1.6/(1+exp(-0.072*(V-5)))
beta_ml(V::Float64) = 0.02*(V-1.31)/(exp((V-1.31)/5.36)-1)
ml_inf(V::Float64) = alpha_ml(V) / (alpha_ml(V) + beta_ml(V))
tau_ml(V::Float64) = 1 / (alpha_ml(V) + beta_ml(V))

# Activación de potasio de tipo A2 / Activation of A2-type Potassium Current
mA2_inf(V::Float64) = 1/(1+exp(-(V+36)/20))
tau_mA2(V::Float64) = 0.37 + 1/(exp((V+35.82)/19.697)+exp((V+79.69)/-12.7))

# Inactivación de potasio de tipo A2 / Inactivation of A2-type Potassium Current
hA2_inf(V::Float64) = 1/(1+exp((V+78)/6))
function tau_hA2(V::Float64)
    if V < -73
        tau_hA2 = 1/(exp((V+46.05)/5)+exp((V+238.4)/-37.45))
    else
        tau_hA2 = 60
    end
    return tau_hA2
end

# Ecuación GHK de la corriente de calcio IT / GHK Equation for T-type Calcium Current
function ITGHK(V::Float64, Ca::Float64)
    z = 2
    F = 96485.3399e-3
    R = 8.314
    T = 310
    Cao = 2e-3

    G = ((z^2*F^2*V)/(R*T))*((Ca-(Cao*exp(-(z*F*V)/(R*T))))/(1-exp(-(z*F*V)/(R*T))))
    return G
end


# Ecuación de la liberación presináptica de transmisor (Destexhe, 1998)
# Presynaptic Neurotransmitter Release Equation (Destexhe, 1998)
function Trans(Vpre::Float64)
    Kp = 5#mV
    Vp = 2#mV 
    Tmax = 0.5#mM

    T = (Tmax)/(1+exp(-(Vpre-Vp)/(Kp) ))
    return T
end


## Función para calcular el voltaje membranal de la interneurona, toma las corrientes intrínsecas (McCormick & Huguenard, 1992) y las sinápticas (Destexhe, 1998)
## Function to compute the membrane voltage of the interneuron, considering intrinsic currents (McCormick & Huguenard, 1992) and synaptic currents (Destexhe, 1998)
function dVint(V::Float64, m::Float64, h::Float64, mp::Float64, mK2::Float64, hK2a::Float64, hK2b::Float64, mc::Float64, mt::Float64, ht::Float64, ml::Float64, mA::Float64, hA::Float64, mA2::Float64, hA2::Float64, Iappi::Float64, Istepi::Float64, gTi::Float64, mAMPA::Float64)
    
(dt)*(1/C)*(-gNai*m^3*h*(V-VNa) - gNapi*mp*(V-VNa) - gK2ai*mK2*hK2a*(V-VK) - gK2bi*mK2*hK2b*(V-VK) - gci*mc*(V-VK) - gTi*mt^2*ht*(V-VCa) - gLi*ml^2*(V-VCa) - gAi*mA^4*hA*(V-VK) - gA2i*mA2^4*hA2*(V-VK) - gNaleaki*(V-VNa) - gKleaki*(V-VK) -gsynAMPA*mAMPA*(V-EsynAMPA)  + Iappi + Istepi)
    
end

## Función para calcular el voltaje membranal de la neurona piramidal, toma las corrientes intrínsecas (McCormick & Huguenard, 1992) y las sinápticas (Destexhe, 1998)
## Function to compute the membrane voltage of the pyramidal neuron, considering intrinsic currents (McCormick & Huguenard, 1992) and synaptic currents (Destexhe, 1998)
function dVpyr(V::Float64, m::Float64, h::Float64, mp::Float64, mK2::Float64, hK2a::Float64, hK2b::Float64, mc::Float64, mt::Float64, ht::Float64, ml::Float64, mA::Float64, hA::Float64, mA2::Float64, hA2::Float64, Iapp::Float64, Istep::Float64, gT::Float64, mGABAA::Float64, rGABAB::Float64, sGABAB::Float64)

(dt)*(1/C)*(-gNa*m^3*h*(V-VNa) - gNap*mp*(V-VNa) - gK2a*mK2*hK2a*(V-VK) - gK2b*mK2*hK2b*(V-VK) - gc*mc*(V-VK) - gT*mt^2*ht*(V-VCa) - gL*ml^2*(V-VCa) - gA*mA^4*hA*(V-VK) - gA2*mA2^4*hA2*(V-VK) - gNaleak*(V-VNa) - gKleak*(V-VK) - gsynGABAA*mGABAA*(V-EsynGABAA) - gsynGABAB*(sGABAB^4)/(sGABAB^4+100)*(V-EKb) + Iapp + Istep)
end


#Ecuaciones diferenciales de las corrientes intrínsecas(McCormick & Huguenard, 1992)
# Differential equations for intrinsic currents (McCormick & Huguenard, 1992)
dm(V::Float64, m::Float64)                 = (dt)* ((1/tau_m(V)) * (m_inf(V) - m))
dh(V::Float64, h::Float64)                 = (dt)* ((1/tau_h(V)) * (h_inf(V) - h))
dmp(V::Float64,mp::Float64)                = (dt)* ((1/tau_mp(V)) * (mp_inf(V) - mp))
dmA(V::Float64, mA::Float64)               = (dt)* ((1/tau_mA(V)) * (mA_inf(V) - mA))
dhA(V::Float64, hA::Float64)               = (dt)* ((1/tau_hA(V)) * (hA_inf(V) - hA))
dmt(V::Float64, mt::Float64)               = (dt)* ((1/tau_mt(V)) * (mt_inf(V) - mt))
dht(V::Float64, ht::Float64)               = (dt)* ((1/tau_ht(V)) * (ht_inf(V) - ht))
dmK2(V::Float64, mK2::Float64)             = (dt)* ((1/tau_mK2(V)) * (mK2_inf(V) - mK2))
dhK2a(V::Float64, hK2a::Float64)           = (dt)* ((1/tau_hK2a(V)) * (hK2a_inf(V) - hK2a))
dhK2b(V::Float64, hK2b::Float64)           = (dt)* ((1/tau_hK2b(V)) * (hK2b_inf(V) - hK2b))
dmc(V::Float64, mc::Float64, CaL::Float64) = (dt)* ((1/tau_mc(V,CaL)) * (mc_inf(V,CaL) - mc))
dml(V::Float64, ml::Float64)               = (dt)* ((1/tau_ml(V)) * (ml_inf(V) - ml))
dCaT(V::Float64, mt::Float64,ht::Float64, CaT::Float64, gT::Float64) = (dt) * (((-5.18e-3*(gT*mt^2*ht*(V-VCa)))/(0.1*29000)) - CaT)  
dCaL(V::Float64, ml::Float64, CaL::Float64) = (dt) * (((-5.18e-3*(gL*ml^2*(V-VCa)))/(0.1*29000)) - CaL); #??
dmA2(V::Float64, mA2::Float64)             = (dt)* ((1/tau_mA2(V)) * (mA2_inf(V) - mA2))
dhA2(V::Float64, hA2::Float64)             = (dt)* ((1/tau_hA2(V)) * (hA2_inf(V) - hA2))


#Ecuaciones diferenciales de las corrientes sinápticas (Destexhe, 1998)
# Differential equations for synaptic currents (Destexhe, 1998)
dmAMPA(V::Float64,mAMPA::Float64) = (dt)*((0.94)*(Trans(V))*(1-mAMPA)-(180*mAMPA))
dmGABAA(V::Float64,mGABAA::Float64) = (dt)*((20)*(Trans(V))*(1-mGABAA)-(160*mGABAA))
dr(V::Float64,rGABAB::Float64) = (dt)*((9)*(Trans(V))*(1-rGABAB)-(1.2*rGABAB))
ds(rGABAB::Float64,sGABAB::Float64) = (dt)*((20*rGABAB)-(17*sGABAB))


#Función en la que se simula la dinámica de la neurona piramidal y la interneurona interconectadas: dVpyr y dVint, respectivamente.
# Function that simulates the dynamics of interconnected pyramidal and interneuron neurons: dVpyr and dVint, respectively.
function simulateTC(Iapp::Float64, Tstepinit::Int64, Tstepfinal::Int64, Istep::Float64, gT::Float64, Iappi::Float64, Istepi::Float64, gTi::Float64)
    # Valores iniciales de la piramidal
    # Initial values for the pyramidal neuron
    V::Float64 = -65.0
    Vprev::Float64 = V
    CaT::Float64 = 50.e-9
    CaL::Float64 = 50.e-9
    m::Float64 = m_inf(V)
    h::Float64 = h_inf(V)
    mp::Float64 = mp_inf(V)
    mA::Float64 = mA_inf(V)
    hA::Float64 = hA_inf(V)
    mt::Float64 = mt_inf(V)
    ht::Float64 = ht_inf(V)
    mK2::Float64 = mK2_inf(V)
    hK2a::Float64 = hK2a_inf(V)
    hK2b::Float64 = hK2b_inf(V)
    mc::Float64 = mc_inf(V,CaL)
    ml::Float64 = ml_inf(V)
    mA2::Float64 = mA2_inf(V)
    hA2::Float64 = hA2_inf(V)
    Vint=zeros(Tdt)
    Vpyr=zeros(Tdt)
    Vpyr[1]=V
    
    #Agregados (Added)
    #---------------------------------------------------------------------------------------------------
    mGABAA_vec = zeros(Tdt)
    rGABAB_vec =zeros(Tdt)
    sGABAB_vec = zeros(Tdt)
    mAMPA_vec = zeros(Tdt)

    #---------------------------------------------------------------------------------------------------
    # Valores iniciales de la interneurona
    # Initial values for the interneuron
    Vi::Float64 = -90.0
    Vprevi::Float64 = Vi
    CaTi::Float64 = 50.e-9
    CaLi::Float64 = 50.e-9
    mi::Float64 = m_inf(Vi)
    hi::Float64 = h_inf(Vi)
    mpi::Float64 = mp_inf(Vi)
    mAi::Float64 = mA_inf(Vi)
    hAi::Float64 = hA_inf(Vi)
    mti::Float64 = mt_inf(Vi)
    hti::Float64 = ht_inf(Vi)
    mK2i::Float64 = mK2_inf(Vi)
    hK2ai::Float64 = hK2a_inf(Vi)
    hK2bi::Float64 = hK2b_inf(Vi)
    mci::Float64 = mc_inf(Vi,CaLi)
    mli::Float64 = ml_inf(Vi)
    mA2i::Float64 = mA2_inf(Vi)
    hA2i::Float64 = hA2_inf(Vi)
    Vint[1]=Vi

    Tstart::Int64 = convert(Int64, round(Tstepinit/dt))
    Tstop::Int64 = convert(Int64, round(Tstepfinal/dt))
    
    #Valores iniciales Agregados (Added initial values)
    mGABAA = 1e-12
    rGABAB = 1e-12    
    sGABAB = 1e-12
    mAMPA = 1e-12

    for z= 2:Tdt
        if z>=Tstart && z<=Tstop
            Iappstep = Istep
        else
            Iappstep = 0.
        end

        #PIRAMIDAL
        # PYRAMIDAL NEURON
        V +=dVpyr(V, m, h, mp, mK2, hK2a, hK2b, mc, mt, ht, ml, mA, hA, mA2, hA2, Iapp, Iappstep, gT, mGABAA, rGABAB, sGABAB)
        m +=dm(Vprev,m)
        h +=dh(Vprev,h)
        mp +=dmp(Vprev,mp)
        mA +=dmA(Vprev,mA)
        hA +=dhA(Vprev,hA)
        mt +=dmt(Vprev,mt)
        ht +=dht(Vprev,ht)
        mK2 +=dmK2(Vprev,mK2)
        hK2a +=dhK2a(Vprev,hK2a)
        hK2b +=dhK2b(Vprev,hK2b)
        mc +=dmc(Vprev,mc, CaL)
        ml +=dml(Vprev,ml)
        CaT += dCaT(Vprev,mt,ht,CaT, gT)
        CaL += dCaL(Vprev,ml,CaL)
        mA2 += dmA2(Vprev,mA2)
        hA2 += dhA2(Vprev,hA2)
        
        #Agregados (Added)
        mGABAA += dmGABAA(Vprevi, mGABAA)
        rGABAB += dr(Vprevi,rGABAB)/100
        sGABAB += ds(rGABAB,sGABAB)/100
        
       
        #---------------------------------------------------------------------------    
        
        #INTERNEURONA
        # INTERNEURON
        Vi +=dVint(Vi, mi, hi, mpi, mK2i, hK2ai, hK2bi, mci, mti, hti, mli, mAi, hAi, mA2i, hA2i, Iappi, Iappstep, gTi, mAMPA)
        mci +=dmc(Vprevi,mci, CaLi)
        CaTi += dCaT(Vprevi,mti,hti,CaTi,gTi)
        CaLi += dCaL(Vprevi,mli,CaLi)
        mi +=dm(Vprevi,mi)
        hi +=dh(Vprevi,hi)
        mpi +=dmp(Vprevi,mpi)
        mAi +=dmA(Vprevi,mAi)
        hAi +=dhA(Vprevi,hAi)
        mti +=dmt(Vprevi,mti)
        hti +=dht(Vprevi,hti)
        mK2i +=dmK2(Vprevi,mK2i)
        hK2ai +=dhK2a(Vprevi,hK2ai)
        hK2bi +=dhK2b(Vprevi,hK2bi)
        mli +=dml(Vprevi,mli)
        mA2i += dmA2(Vprevi,mA2i)
        hA2i += dhA2(Vprevi,hA2i)
        
        #Agregado (added)
        mAMPA += dmAMPA(Vprev, mAMPA)
        
        Vprev = copy(V)
        Vprevi = copy(Vi)
        
        Vpyr[z] = copy(V)
        Vint[z] = copy(Vi)
        mGABAA_vec[z] = copy(mGABAA)
        rGABAB_vec[z] = copy(rGABAB)
        sGABAB_vec[z] = copy(sGABAB)
        mAMPA_vec[z] = copy(mAMPA)

    end
    return [Vpyr, Vint, mGABAA_vec, sGABAB_vec, rGABAB_vec, mAMPA_vec]
end
