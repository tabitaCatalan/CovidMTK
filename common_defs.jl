dt = 1. # Intervalo de sampleo de observaciones 
T = 500. # Se resolverá el problema en el intervalo [0,T]
if is_synthetic
    T = 230.
end
Nmediciones = Int(T/dt) # número de mediciones 
ts = 0.:dt:(T-dt) # grilla de tiempos 

t0date = Date(2020,3,30)
tsdate = t0date:Day(1):(t0date + Day(length(ts) - 1))