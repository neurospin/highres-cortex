# coding: utf-8
from soma import aims
from soma import aimsalgo

gradnx = aims.read("heat_gradnx.nii")
gradny = aims.read("heat_gradny.nii")
gradnz = aims.read("heat_gradnz.nii")

gradient = aimsalgo.AimsGradient_FLOAT()
div_gradn = (gradient.X(gradnx) +
             gradient.Y(gradny) +
             gradient.Z(gradnz))
aims.write(div_gradn, "./heat_div_gradn.nii.gz")
