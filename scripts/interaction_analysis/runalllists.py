import os
import sys, argparse

print("Running HH-GTEx-CAN")
os.system('python interactions.py --subset hedgehog --dataset GTEx --set CAN')

print("Running HH-GTEx-NONCAN")
os.system('python interactions.py --subset hedgehog --dataset GTEx --set NONCAN')

print("Running HH-GTEx-RAN")
os.system('python interactions.py --subset hedgehog --dataset GTEx --set RANDOM')

print("Running HH-TCGA-CAN")
os.system('python interactions.py --subset hedgehog --dataset TCGA --set CAN')

print("Running HH-TCGA-NONCAN")
os.system('python interactions.py --subset hedgehog --dataset TCGA --set NONCAN')

print("Running HH-TCGA-RAN")
os.system('python interactions.py --subset hedgehog --dataset TCGA --set RANDOM')

print("Running NOTCH-GTEx-CAN")
os.system('python interactions.py --subset notch --dataset GTEx --set CAN')

print("Running NOTCH-GTEx-RAN")
os.system('python interactions.py --subset notch --dataset GTEx --set RANDOM')

print("Running NOTCH-GTEx-NONCAN")
os.system('python interactions.py --subset notch --dataset GTEx --set NONCAN')

print("Running NOTCH-TCGA-CAN")
os.system('python interactions.py --subset notch --dataset TCGA --set CAN')

print("Running NOTCH-TCGA-NONCAN")
os.system('python interactions.py --subset notch --dataset TCGA --set NONCAN')

print("Running NOTCH-TCGA-RAN")
os.system('python interactions.py --subset notch --dataset TCGA --set RANDOM')


print("Running ANGIOGENESIS-GTEx-CAN")
os.system('python interactions.py --subset angiogenesis --dataset GTEx --set CAN')

print("Running ANGIOGENESIS-GTEx-NONCAN")
os.system('python interactions.py --subset angiogenesis --dataset GTEx --set NONCAN')

print("Running ANGIOGENESIS-GTEx-RAN")
os.system('python interactions.py --subset angiogenesis --dataset GTEx --set RANDOM')

print("Running ANGIOGENESIS-TCGA-CAN")
os.system('python interactions.py --subset angiogenesis --dataset TCGA --set CAN')

print("Running ANGIOGENESIS-TCGA-NONCAN")
os.system('python interactions.py --subset angiogenesis --dataset TCGA --set NONCAN')

print("Running ANGIOGENESIS-TCGA-RAN")
os.system('python interactions.py --subset angiogenesis --dataset TCGA --set RANDOM')

print("Running MYCV2-GTEx-CAN")
os.system('python interactions.py --subset mycv2 --dataset GTEx --set CAN')

print("Running MYCV2-GTEx-NONCAN")
os.system('python interactions.py --subset mycv2 --dataset GTEx --set NONCAN')

print("Running MYCV2-GTEx-RAN")
os.system('python interactions.py --subset mycv2 --dataset GTEx --set RANDOM')

print("Running MYCV2-TCGA-CAN")
os.system('python interactions.py --subset mycv2 --dataset TCGA --set CAN')

print("Running MYCV2-TCGA-NONCAN")
os.system('python interactions.py --subset mycv2 --dataset TCGA --set NONCAN')

print("Running MYCV2-TCGA-RAN")
os.system('python interactions.py --subset mycv2 --dataset TCGA --set RANDOM')
