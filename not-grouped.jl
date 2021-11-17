comunas2 = [
    13123, # Providencia 
    #13110, # La Florida 
    #13101, # Santiago 
    #13201, # Puente Alto 
    #13119, # Maipú 
    #13106, # Estación Central 
    #13126, # Quinta Normal 
    13401, # San Bernardo 
    #13105, # El Bosque 
    #13112, # La Pintana 
    #13116, # Lo Espejo 
]
#comunas2 = rand(comunas, number_of_selected_municipalities)
    
Pt = makePmatrix(initial_mob, make_mob_in_time(comunas2, datamap)); 