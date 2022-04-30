function ambient_to_properties(f=1000.0,ps=101325.0,t=20.0,H=50.0)

    # [rho,c0,cf,CpCv,eta,alfa,Cp,Cv,landa,beta]=ambient_to_properties(ps,t,H,f)
    #
    # Calculation of the physical properties of air ambient conditions.
    # Formulas and notation after Rasmussen [1] (annex A).
    #
    # Input:
    #  -ps:     static pressure in Pa. Standard: 101325 Pa
    #  -t:      temperature in º Celsius. Standard: 20 Celsius
    #  -H:      humidity in #. Standard: 50#
    #  -f:      frequency in Hz
    #
    # Output:
    #  -rho:    Density of air, kg/m3
    #  -c0:     Zero-frequency speed of sound in air, m/s
    #  -cf:     Speed of sound at actual frequency, m/s
    #  -CpCv:   Ratio of specific heats
    #  -eta:    Viscosity of air, N*s/m2 = Pa*s. It is called "mu" elsewhere. 
    #  -alfa:   Air attenuation coefficient in Np/m
    #  -Cp:     Specific heat capacity at constant pressure, J/(kg*K)
    #  -Cv:     Specific heat capacity at constant volume, J/(kg*K)
    #  -landa:  Thermal conductivity, J/(s*m*K)=W/(m*K)
    #  -beta    =rho*(Cp-Cv); thermal expansion coefficient, Pa/K
    #
    #
    # [1] Knud Rasmussen. "Calculation methods for the physical properties of air
    #     used in the calibration of microphones". Report PL-11b, Department of
    #     Acoustic Technology, Technical University of Denmark, 1997.
    
    # Matlab version by V. Cutanda 05-2006
    # Added outputs, 02-2007, V. Cutanda
    
    
    T0=273.15;          # 0 ºC in K
    T20=293.15;         # 20 ºC in K
    psr=101325;         # Reference static pressure
    T=T0+t;             # Temperature in K
    
    
    # Matrix of coefficients in table A.1 in [1], columns correspond to:
    #   1: Saturation water vapor pressure psv
    #   2: Enhancement factor f(ps,t)
    #   3: Compressibility factor Z
    #   4: Zero-frequency speed of sound c0
    #   5: Ratio of specific heats ?
    #   6: Viscosity ?
    #   7: Thermal conductivity landa
    #   8: Specific heat capacity at constant pressure Cp
    aa= [1.2378847e-5 1.00062 1.58123e-6 331.5024 1.400822 84.986 60.054 0.251625;   
         -1.9121316e-2 3.14e-8 -2.9331e-8 0.603055 -1.75e-5 7 1.846 -9.2525e-5;    
         33.93711047 5.6e-7 1.1043e-10 -0.000528 -1.73e-7 113.157 2.06e-6 2.1334e-7; 
         -6.3431645e3 NaN 5.707e-6 51.471935 -0.0873629 -1 40 -1.0043e-10; 
         NaN NaN -2.051e-8 0.1495874 -0.0001665 -3.7501e-3 -1.775e-4 0.12477; 
         NaN NaN 1.9898e-4 -0.000782 -3.26e-6 -100.015 NaN -2.283e-5; 
         NaN NaN -2.376e-6 -1.82e-7 2.047e-8 NaN NaN 1.267e-7; 
         NaN NaN 1.83e-11 3.73e-8 -1.26e-10 NaN NaN 0.01116; 
         NaN NaN -0.765e-8 -2.93e-10 5.939e-14 NaN NaN 4.61e-6; 
         NaN NaN NaN -85.20931 -0.1199717 NaN NaN 1.74e-8; 
         NaN NaN NaN -0.228525 -0.0008693 NaN NaN NaN; 
         NaN NaN NaN 5.91e-5 1.979e-6 NaN NaN NaN; 
         NaN NaN NaN -2.835149 -0.01104 NaN NaN NaN; 
         NaN NaN NaN -2.15e-13 -3.478e-16 NaN NaN NaN; 
         NaN NaN NaN 29.179762 0.0450616 NaN NaN NaN; 
         NaN NaN NaN 0.000486 1.82e-6 NaN NaN NaN];
    
      
    # Saturation water vapor pressure:
    # psv = exp(aa[1,1]*T^2 + aa[2,1]*T + aa[3,1) + aa[4,1)/T); # Formula Giacomo
    psv = psr * 10^(4.6151 - 6.8346*((T0+0.01)/T)^1.261);     # Formula ISO 9613-1:1993
    
    # Enhancement factor:
    fpst = aa[1,2] + aa[2,2]*ps + aa[3,2]*t^2;
    
    # Mole fraction of water vapor in air:
    xw = H*psv*fpst/(100*ps);
    
    # Compressibility factor:
    Z = 1 -  ps/T*(aa[1,3] + aa[2,3]*t + aa[3,3]*t^2 + (aa[4,3]+aa[5,3])*xw + (aa[6,3]+aa[7,3])*xw^2) + 
         (ps/T)^2*(aa[8,3] + aa[9,3]*xw^2);
      
    # Mole fraction of carbon-dioxide in air:
    xc=0.0004;
    
    # Density of air:
    rho = (3.48349 + 1.44*(xc-0.0004)) * 1e-3*ps/(Z*T) * (1-0.3780*xw);
    
    # Zero-frequency speed of sound in air:
    c0 =  aa[1,4] + aa[2,4]*t + aa[3,4]*t^2 + (aa[4,4] + aa[5,4]*t + aa[6,4]*t^2)*xw + 
          (aa[7,4] + aa[8,4]*t + aa[9,4]*t^2)*ps + (aa[10,4] + aa[11,4]*t + aa[12,4]*t^2)*xc + 
          aa[13,4]*xw^2 + aa[14,4]*ps^2 + aa[15,4]*xc^2 + aa[16,4]*xw*ps*xc;
    
    # Ratio of specific heats:
    CpCv = aa[1,5] + aa[2,5]*t + aa[3,5]*t^2 + (aa[4,5] + aa[5,5]*t + aa[6,5]*t^2)*xw + 
          (aa[7,5] + aa[8,5]*t + aa[9,5]*t^2)*ps + (aa[10,5] + aa[11,5]*t + aa[12,5]*t^2)*xc + 
           aa[13,5]*xw^2 + aa[14,5]*ps^2 + aa[15,5]*xc^2 + aa[16,5]*xw*ps*xc;
       
    # Viscosity of air:
    eta = (aa[1,6] + aa[2,6]*T + (aa[3,6]+aa[4,6]*T)*xw + aa[5,6]*T^2 + aa[6,6]*xw^2) * 1e-8;
    
    # Thermal conductivity:
    landa = (aa[1,7] + aa[2,7]*T + aa[3,7]*T^2 + (aa[4,6] + aa[5,6]*T)*xw) * 1e-8;
    
    # Specific heat capacity at constant pressure
    Cp = aa[1,8] + aa[2,8]*T + aa[3,8]*T^2 + aa[4,8]*T^3 + (aa[5,8] + aa[6,8]*T + aa[7,8]*T^2)*xw +
         (aa[8,8] + aa[9,8]*T + aa[10,8]*T^2)*xw^2;
    
    # Diffusivity of air;
    alfa_t =  landa/(rho*Cp);
    
    # Relaxation frequency of oxigen:
    frO = ps/psr* (24 + 4.04e6*xw*(0.2+1e3*xw)/(3.91+1e3*xw));
    
    # Relaxation frequency of nitrogen:
    frN = ps/psr* (T/T20)^(-1/2) * (9 + 28.0e3*xw*exp(-4.170*((T/T20)^(-1/3)-1)));
    
    # Attenuation coefficient of relaxation in oxygen in Np/m:
    alfa_vO = 0.01275*f^2*exp(-2239.1/T)/(frO+f^2/frO)*(T20/T)^(5/2);
    
    # Attenuation coefficient of relaxation in nitrogen in Np/m:
    alfa_vN = 0.1068*f^2*exp(-3352.0/T)/(frN+f^2/frN)*(T20/T)^(5/2);
    
    # Speed of sound at actual frequency:
    cf = (1/c0 - alfa_vO/(2*pi*frO) - alfa_vN/(2*pi*frN))^-1;
    
    # Air attenuation coefficient in Np/m:
    alfa = 18.4e-12 * f^2 * (ps/psr)^(-1) * (T/T20)^(1/2) + alfa_vO + alfa_vN;
    
    
    # Conversion to SI units:
    Cp=Cp*4.1868*1000;          # from cal/(g*K) to J/(kg*K)
    landa=landa*4.1868*1000;    # from Kcal/(s*m*K) to J/(s*m*K) . The PL11 report wrongly (?) says "cal/g"
    
    # Derived parameters:
    Cv   = Cp/CpCv;
    beta = rho*(Cp-Cv);

    return rho,c0,cf,CpCv,eta,alfa,Cp,Cv,landa,beta
end
    


function betag(sigma,f)

    # betaDB=betag(sigma,f);
    #
    # Calculates the admittance by Delany and Bazley formula
    #
    # Input:
    #    -sigma: flow resistance in Nsm-4. It may be a vector.
    #    -f:     frequency. It may be a vector.
    #
    # Output variables:
    #    -betaDB:  normalised admittance. As many rows as the length of sigma.
    #              As many columns as the length of f.
    
    # Reference: Propagation of Sound in Porous Media: Modelling Sound Absorbing Materials,
    # Second Edition, Jean F. Allard, Noureddine Atalla. 2009, John Wiley.
    # Section 2.5.3
    
    # Modified with new expression, Vicente Cutanda Henríquez 4-2011.
    
    
    # Ambient conditions
    pa = 101325;         # Static pressure (Pa)
    t = 20;              # Temperature (ºC)
    Hr = 50;             # Relative humidity (#)
    rho,c,cf,CpCv,nu,alfa=ambient_to_properties(pa,t,Hr,1000); 
    
    
    betaDB = zeros(length(sigma),length(f));
    
    for ff=1:length(f)
        for rr=1:length(sigma)
            if isinf(sigma[rr])
                betaDB[rr,ff]=0.0;
            elseif isnan(sigma[rr])
                betaDB[rr,ff]=NaN;
            else
                #betaDB(rr,ff)=rho*c./(1+9.08*(1000*f(ff)./sigma(rr)).^(-0.75)+11.9i.*(1000*f(ff)./sigma(rr)).^(-0.73));
                X=rho*f[ff]./sigma[rr];
                betaDB[rr,ff]=1.0 ./ (1+0.0571*(X.^-0.754)+j*0.087*(X.^-0.732)); # normalized admittance, exp(-jwt) convention
                # betaDB(rr,ff)=betaDB(rr,ff)/(rho*c); # de-normalize
            end
        end
    end
    return betaDB
end
