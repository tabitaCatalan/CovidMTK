dt = 1.#0.1 # Intervalo de sampleo de observaciones 
T = 400. # Se resolverá el problema en el intervalo [0,T]
Nmediciones = Int(T/dt) # número de mediciones 
ts = 0.:dt:(T-dt) # grilla de tiempos 

t0date = Date(2020,3,30)
tsdate = t0date:Day(1):(t0date + Day(length(ts) - 1))