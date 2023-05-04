# NUOPC_FVCOM

src -- FVCOM_source code 5.0

nuopc_fvcom_cap -- NUOPC FVCOM Cap source code

wmesmfmd.ftn -- Added code to get significant wave height, average wave length and
                average wave direction for export to coupler
		
module_EARTH_GRID_COMP.F90 -- Added code to allow NEMS to implenment NUOPC FVCOM Cap
                and added fields from WW3 to FVCOM :
                'significant_wave_height'
                'average_wave_length'
                'average_wave_direction'	
