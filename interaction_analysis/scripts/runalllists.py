import os
import sys, argparse

print("Running HH-GTEx-CAN")
os.system('python interactions.py --subset hedgehog --dataset GTEx --set CAN')

print("Running HH-GTEx-NONCAN")
os.system('python interactions.py --subset hedgehog --dataset GTEx --set NONCAN')

print("Running HH-TCGA-CAN")
os.system('python interactions.py --subset hedgehog --dataset TCGA --set CAN')

print("Running HH-TCGA-NONCAN")
os.system('python interactions.py --subset hedgehog --dataset TCGA --set NONCAN')

#TODO RANDOMS 5 ITERS

print("Running NOTCH-GTEx-CAN")
os.system('python interactions.py --subset notch --dataset GTEx --set CAN')

print("Running NOTCH-GTEx-NONCAN")
os.system('python interactions.py --subset notch --dataset GTEx --set NONCAN')

print("Running NOTCH-TCGA-CAN")
os.system('python interactions.py --subset notch --dataset TCGA --set CAN')

print("Running NOTCH-TCGA-NONCAN")
os.system('python interactions.py --subset notch --dataset TCGA --set NONCAN')

#TODO RANDOMS 5 ITERS

print("Running ANGIOGENESIS-GTEx-CAN")
os.system('python interactions.py --subset angiogenesis --dataset GTEx --set CAN')

print("Running NOTCH-ANGIOGENESIS-NONCAN")
os.system('python interactions.py --subset angiogenesis --dataset GTEx --set NONCAN')

print("Running NOTCH-ANGIOGENESIS-CAN")
os.system('python interactions.py --subset angiogenesis --dataset TCGA --set CAN')

print("Running NOTCH-ANGIOGENESIS-NONCAN")
os.system('python interactions.py --subset angiogenesis --dataset TCGA --set NONCAN')

#TODO RANDOMS 5 ITERS

print("Running MYCV2-GTEx-CAN")
os.system('python interactions.py --subset mycv2 --dataset GTEx --set CAN')

print("Running MYCV2-ANGIOGENESIS-NONCAN")
os.system('python interactions.py --subset mycv2 --dataset GTEx --set NONCAN')

print("Running MYCV2-ANGIOGENESIS-CAN")
os.system('python interactions.py --subset mycv2 --dataset TCGA --set CAN')

print("Running MYCV2-ANGIOGENESIS-NONCAN")
os.system('python interactions.py --subset mycv2 --dataset TCGA --set NONCAN')

#TODO RANDOMS 5 ITERS
