
#cmd.select( 'flexchain', 'resi 22-30' )
#cmd.hide( 'cartoon', 'flexchain' )
#cmd.select( 'a_domain', 'resi 31-140' )
#cmd.color( 'blue', 'a_domain' )
#cmd.select( 'b_domain', 'resi 141-240' )
#cmd.color( 'green', 'b_domain' )
#cmd.select( 'bprime_domain', 'resi 241-359' )
#cmd.color( 'yellow', 'bprime_domain' )
#cmd.select( 'x_linker', 'resi 360-374' )
#cmd.color( 'black', 'x_linker' )
#cmd.select( 'aprime_domain', 'resi 375-484' )
#cmd.color( 'orange', 'aprime_domain' )
#cmd.select( 'c_domain', 'resi 485-504' )
#cmd.color( 'red', 'c_domain' )
#cmd.select( 'residue61', 'resi 61 and name ca')
#cmd.color( 'yellow', 'residue61')
#cmd.show( 'sphere', 'residue61' )
#cmd.select( 'residue406', 'resi 406 and name ca')
#cmd.color( 'yellow', 'residue406')
#cmd.show( 'sphere', 'residue406' )


#cmd.set(full_screen='on')

### cut below here and paste into script ###
# cmd.set_view('0.979684353,   -0.145519406,   -0.137993217,\
#      0.135169536,   -0.029161714,    0.990392923,\
#     -0.148145691,   -0.988925874,   -0.008899692,\
#      0.000000000,    0.000000000, -348.124450684,\
#    -14.108730316,  -18.215091705,   66.387222290,\
#    274.463958740,  421.784942627,  -20.000000000' )

### cut above here and paste into script ###
cmd.hide()
cmd.show('cartoon')

cmd.bg_color('grey')
cmd.movie.add_state_sweep(2,1,'start=1')
cmd.movie.produce(filename,mode='ray',quality=100,quiet=1,preserve=1)
cmd.save(folder+'cartoon.pse')

print ('CARTOON rendering finished')

cmd.quit()
