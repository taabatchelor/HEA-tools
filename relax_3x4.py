'''relaxes the atomic structures and saves them to a database'''

# load python modules
from gpaw import GPAW, PW, FermiDirac
from ase.optimize import QuasiNewton
from ase.db import connect

# define the database containing unrelaxed slabs
db = connect('../1_make_surfaces/3x4.db')

# define the database to be filled with relaxed structures
db_relaxed = connect('test_trial.db')

# loop through all 1000 unrelaxed structures
for idx in range(1000):
	
	## slab
	# get slab
	slab = db.get_atoms(relaxed=0, type='slab', slabId=idx)
	
	# define slab calculator
	calc = GPAW(mode = PW(400),					# PW energy cutoff is 400 eV
				basis = 'dzp',					# double zeta potential basis
				xc="RPBE",						# exchange-correlation functional
				eigensolver = 'rmm-diis',		# rmm-diis eigensolver
				kpts=(4, 4, 1),					# k-points in (x, y, z)-directions
				occupations = FermiDirac(0.1),	# Fermi Dirac occupasions set to 0.1
				txt="slab_%d.txt"%idx)			# output text file

	# attach calculator to the atoms object
	slab.set_calculator(calc)

	# relax slab to a maximum force of 0.1 eV/A
	dyn = QuasiNewton(slab, trajectory="slab_%d.traj"%idx)
	dyn.run(fmax=0.1)

	# write relaxed slab to the database
	db_relaxed.write(slab, relaxed=1, type='slab', slabId=idx)
	
	
	## OH
	# get slab with OH adsorbed
	slabOH = db.get_atoms(relaxed=0, type='OH', slabId=idx)
	
	# define slab + OH calculator
	calc = GPAW(mode = PW(400),					# PW energy cutoff is 400 eV
				basis = 'dzp',					# double zeta potential basis
				xc="RPBE",						# exchange-correlation functional
				eigensolver = 'rmm-diis',		# rmm-diis eigensolver
				kpts=(4, 4, 1),					# k-points in (x, y, z)-directions
				occupations = FermiDirac(0.1),	# Fermi Dirac occupasions set to 0.1
				txt="OH_%d.txt"%idx)			# output text file
	
	# attach calculator to the atoms object
	slabOH.set_calculator(calc)
	
	# relax slab with OH to a maximum force of 0.1 eV/A
	dyn = QuasiNewton(slabOH, trajectory="OH_%d.traj"%idx)
	dyn.run(fmax=0.1)
	
	# write relaxed slab with OH to the database
	db_relaxed.write(slabOH, relaxed=1, type='OH', slabId=idx)
	
	
	## O
	# get slab with O adsorbed
	slabO = db.get_atoms(relaxed=0, type='O', slabId=idx)
	
	# define slab + O calculator
	calc = GPAW(mode = PW(400),					# PW energy cutoff is 400 eV
				basis = 'dzp',					# double zeta potential basis
				xc="RPBE",						# exchange-correlation functional
				eigensolver = 'rmm-diis',		# rmm-diis eigensolver
				kpts=(4, 4, 1),					# k-points in (x, y, z)-directions
				occupations = FermiDirac(0.1),	# Fermi Dirac occupasions set to 0.1
				txt="O_%d.txt"%idx)				# output text file
	
	# relax slab with O to a maximum force of 0.1 eV/A
	dyn = QuasiNewton(slabO, trajectory="O_%d.traj"%idx)
	dyn.run(fmax=0.1)

	# write relaxed slab with OH to the database
	db_relaxed.write(slabO, relaxed=1, type='O', slabId=idx)
