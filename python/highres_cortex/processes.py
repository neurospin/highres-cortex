import subprocess

from capsul.process import Process
from traits.api import File, Float, Int, Undefined

class LaplacianProcess(Process):
    classif = File(Undefined, output=False)
    precision = Float(0.001, output=False)
    typical_cortical_thickness = Float(3, output=False)
    verbosity = Int(1, output=False)

    laplace_field = File(Undefined, output=True)

    def get_commandline(self):
        return ["ylLaplacian",
                "--classif", self.classif,
                "--output", self.laplace_field,
                "--precision", self.precision,
                "--typical-cortical-thickness", self.typical_cortical_thickness,
                "--verbose", self.verbosity]
