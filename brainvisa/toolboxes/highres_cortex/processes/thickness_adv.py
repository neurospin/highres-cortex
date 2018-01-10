from brainvisa.processes import *
from brainvisa.processing import capsul_process

name = "Cortical thickness (advection)"
userLevel = 1

base_class = capsul_process.CapsulProcess
capsul_process = "highres_cortex.capsul.thickness_adv"
