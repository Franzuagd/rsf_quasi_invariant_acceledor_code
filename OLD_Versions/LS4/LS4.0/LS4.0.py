######################## IMPORTS #####################################

import subprocess
import re
import os
#os.environ["OMP_NUM_THREADS"] = "1"
#os.environ["OPENBLAS_NUM_THREADS"] = "1"
#os.environ["MKL_NUM_THREADS"] = "1"
#os.environ["NUMEXPR_NUM_THREADS"] = "1"
import multiprocessing
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import time
import scipy
import random
import pygad
from typing import List
from functools import partial
from scipy.stats import skew
from numpy import linalg as LA
from scipy import linalg
from scipy.interpolate import make_interp_spline
from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)
random.seed(42)
import sympy as sp
import numpy as np
import sympy as sp
from math import comb
import scipy as sc
import scipy.sparse as sps
from collections import defaultdict

np.set_printoptions(
    precision=17,
    suppress=False,
    linewidth=100000,
    threshold=np.inf
)

delta, x, y, px, py = sp.symbols('delta x y px py')
b1, b2, b3, b4, b5 = sp.symbols('b1 b2 b3 b4 b5')
vars = [delta, x, y, px, py]

# Default tolerances and limits. main() can overwrite these values.
nulm2e_tol = 1e-12
invariant_lstsq_tol = 1e-12
invariant_bracket_zero_tol = 1e-25
norm_floor_tol = 1e-300
max_natemit = 71e-12
min_dispersion = 0.0
max_initial_dispersion = 1e-2
max_beta_y = 22
max_alpha_x = 1e-3
max_linear_trace = 1.9
linear_trace_penalty_scale = 1e28


#################################################################
#OLD CODE FUNCTIONS  
##################################################################

########### Read files ##########

def oparing(opafile):
	with open(opafile) as f:
		lines = f.readlines()
#######################Reading and extracting information from OPA file#######################
	lines1=[]
	fc=0
	for i in range(0,len(lines)):
		if len(lines[i].strip())>1:
			if fc==0:
				if lines[i].strip()[len(lines[i].strip())-1]==';':
					lines1.append(lines[i])
				elif lines[i].strip()[len(lines[i].strip())-1]==',':
					lines2=lines[i].strip()
					fc=1
			else:
				if lines[i].strip()[len(lines[i].strip())-1]==',':
					lines2=lines2+lines[i].strip()
				else:
					lines2=lines2+lines[i].strip()
					lines1.append(lines2)
					fc=0
					del(lines2)
	elements=pd.DataFrame()
	variables=[]
	cells=[]
	for i in lines1:
		en=i.strip().lower().split('=',1)
		if en[0].strip()=='energy':
			energy=float(en[1].strip().lower().split(';',1)[0])
		a=i.strip().split(':', 1)
		if len(a)>1:
			b=a[1].strip().split(',',1)
			if b[0].lower()=='drift':
				c=b[1].strip().lower().split('l',1)
				d=c[1].strip().split(',',1)
				e=d[0].strip().split('=',1)
				elements[a[0].strip().lower()]=['drift',str(e[1].strip()),0,0,0,0,0,'0']
			elif b[0].lower()=='bending':
				c=b[1].strip().lower().split('l',1)
				d=c[1].strip().split(',',1)
				e=d[0].strip().split('=',1)
				f=b[1].strip().lower().split('t',1)
				g=f[1].strip().split(',',1)
				h=g[0].strip().split('=',1)
				f1=b[1].strip().lower().split('k',1)
				g1=f1[1].strip().split(',',1)
				h1=g1[0].strip().split('=',1)
				f2=b[1].strip().lower().split('t1',1)
				g2=f2[1].strip().split(',',1)
				h2=g2[0].strip().split('=',1)
				f3=b[1].strip().lower().split('t2',1)
				g3=f3[1].strip().split(',',1)
				h3=g3[0].strip().split('=',1)
				elements[a[0].strip().lower()]=['bending',str(e[1].strip()),str(h[1].strip()),str(h1[1].strip()),str(h2[1].strip()),str(h3[1].strip()),0,'1']
			elif b[0].lower()=='quadrupole':
				c=b[1].strip().lower().split('l',1)
				d=c[1].strip().split(',',1)
				e=d[0].strip().split('=',1)
				f1=b[1].strip().lower().split('k',1)
				g1=f1[1].strip().split(',',1)
				h1=g1[0].strip().split('=',1)
				elements[a[0].strip().lower()]=['quadrupole',str(e[1].strip()),0,str(h1[1].strip()),0,0,0,'2']
			elif b[0].lower()=='sextupole':
				c=b[1].strip().lower().split('l',1)
				d=c[1].strip().split(',',1)
				e=d[0].strip().split('=',1)
				f1=b[1].strip().lower().split('k',1)
				g1=f1[1].strip().split(',',1)
				h1=g1[0].strip().split('=',1)
				elements[a[0].strip().lower()]=['sextupole',str(e[1].strip()),0,0,0,0,str(h1[1].strip()),'3']
			elif b[0].lower()=='multipole':
				c=b[1].strip().lower().split('n',1)
				d=c[1].strip().split(',',1)
				e=d[0].strip().split('=',1)
				f1=b[1].strip().lower().split('k',1)
				g1=f1[1].strip().split(',',1)
				h1=g1[0].strip().split('=',1)
				elements[a[0].strip().lower()]=['multipole',str('1e-8'),0,0,0,0,0,str(h1[1].strip())]
			elif b[0].lower()=='octupole':
				c=b[1].strip().lower().split('n',1)
				d=c[1].strip().split(',',1)
				e=d[0].strip().split('=',1)
				f1=b[1].strip().lower().split('k',1)
				g1=f1[1].strip().split(',',1)
				h1=g1[0].strip().split('=',1)
				elements[a[0].strip().lower()]=['multipole',str('1e-8'),0,0,0,0,0,str(h1[1].strip())]
			else:
				d=a[1].strip().lower().split(',')
				d=[n.strip() for n in d]
				cells.append([a[0].strip().lower().split(','),d])
		else:
			b=a[0].strip().split('=',1)
			if len(b)>1:
				d=b[1].strip().split(';',1)
				variables.append([b[0].strip().lower(),d[0]])
			else:
				b=a[0].strip().split('{',1)
				if len(b)==1 and i!='\n':
					d=a[0].strip().lower().split(',')
					d=[n.strip() for n in d]
					cells.append([d])
	elements.index=['type','length','T0','K','T1','T2','S','N']
	for i in variables:
		if i[0]=='f':
			elements=elements.replace('*'+i[0],'*'+i[1])
		else:
			elements=elements.replace(i[0],i[1])
#		print(i[0],i[1])
#	for j in variables:
#		print(j[0],j[1])
#	for i in elements.columns:
#		print(elements[i])
	for j in variables:
		for i in elements.columns:
#			print(j[0],j[1])
#			print(elements[i])
			elements[i]=elements[i].str.replace('-'+j[0],'-'+str(j[1]))
			elements[i]=elements[i].str.replace('*'+j[0],'*'+str(j[1]))
#			elements[i]=elements[i].str.replace('\*'+j[0],'*'+str(j[1]))
	elements=elements.fillna('2-2')
#	print(elements)
#	elements=map(evaluate_cell,elements)
	elements=evaluate_dataframe(elements.copy())
#	print(elements)
	ring=[]
	for i in range(0,len(cells)):
		for j in range(0,len(cells[len(cells)-1-i][1])):
			a=cells[len(cells)-1-i][1][j].split('*')
			for k in a:
				if k.isdigit():
					[ring.append(a[1].split(';')[0]) for x in range(0,int(k))]
					cells[len(cells)-1-i][1][j]=ring
					ring=[]
	ring2=[]
	ring=[]
	kk=1
	while kk>0:
		for i in reversed(cells):
			kk=0
			for j in i[1]:
				if type(j) is list:
					kk=kk+1
					for k in j:
						ring2.append(k)
				else:
					ring2.append(j.split(';')[0])
			ring.append([i[0],ring2])
			ring2=[]
	for ii in range(0,len(ring)):
		for jj in range(0,len(ring[ii][1])):
			for k in ring:
				if ring[ii][1][jj]==k[0][0]:
					ring[ii][1][jj]=k[1]
	ring3=[]
	ringf=[]
	kk=1
	while kk>0:
		kk=0
		for i in ring:
			for j in i[1]:
				if type(j) is list:
					kk=kk+1
					for k in j:
						ring3.append(k)
				else:
					ring3.append(j.split(';')[0])
			ringf.append([i[0],ring3])
			ring3=[]
		ring=ringf
		ringf=[]
	ring.reverse()
	kj=1
	while kj>0:
		kj=0
		for i in range(0,len(ring)):
			for j in range(0,len(ring[i][1])):
				jj=ring[i][1][j]
				if jj[0]=='-':
					kj=kj+1
					for k in range(0,len(ring)):
						if ring[k][0][0]==jj[1:len(jj)]:
							ring4=[]
							for l in range(0,len(ring[k][1])):
								ring4.append(ring[k][1][len(ring[k][1])-1-l])
							ring[i][1][j]=ring4
				for iii in ring:
					for jjj in iii[1]:
						if type(jjj) is list:
							for kkk in jjj:
								ring3.append(kkk)
						else:
							ring3.append(jjj.split(';')[0])
					ringf.append([iii[0],ring3])
					ring3=[]
				ring=ringf
				ringf=[]
		ring3=[]
		ringf=[]
		kk=1
		while kk>0:
			kk=0
			for i in ring:
				for j in i[1]:
					if type(j) is list:
						kk=kk+1
						for k in j:
							ring3.append(k)
					else:
						ring3.append(j.split(';')[0])
				ringf.append([i[0],ring3])
				ring3=[]
			ring=ringf
			ringf=[]
	return ring,elements,energy

############## General propose ############

def extractelem(lst,elem):
    return [item[elem] for item in lst]

def evaluate_cell(cell):
	try:
		return eval(cell)
	except (SyntaxError, NameError, TypeError):
		return cell

def evaluate_dataframe(df):
	return df.applymap(evaluate_cell)

############# Linear Functions Tranfer Matrices ###############

def driff(L):
	matrix= np.array([[1,L,0,0],[0,1,0,0],[0,0,1,L],[0,0,0,1]])
	return matrix

def q(L,K):
	if K>0:
		matrix= np.array([[math.cos(L*math.sqrt(K)),math.sin(L*math.sqrt(K))/math.sqrt(K),0,0],[-math.sqrt(K)*math.sin(L*math.sqrt(K)),math.cos(L*math.sqrt(K)),0,0],[0,0,math.cosh(L*math.sqrt(K)),math.sinh(L*math.sqrt(K))/math.sqrt(K)],[0,0,math.sqrt(K)*math.sinh(L*math.sqrt(K)),math.cosh(L*math.sqrt(K))]])
	elif K<0:
		K=-K
		matrix= np.array([[math.cosh(L*math.sqrt(K)),math.sinh(L*math.sqrt(K))/math.sqrt(K),0,0],[math.sqrt(K)*math.sinh(L*math.sqrt(K)),math.cosh(L*math.sqrt(K)),0,0],[0,0,math.cos(L*math.sqrt(K)),math.sin(L*math.sqrt(K))/math.sqrt(K)],[0,0,-math.sqrt(K)*math.sin(L*math.sqrt(K)),math.cos(L*math.sqrt(K))]])
	else:
		matrix= np.array([[1,L,0,0],[0,1,0,0],[0,0,1,L],[0,0,0,1]])
	return matrix

def cfsd(L,K,T0,L0):
	h=(np.pi*T0)/(180*L0)
	KH=(K+h**2)**0.5
	if K==0:
		matrix=np.array([[np.cos(L*h),np.sin(L*h)/h,0,0],[-h*np.sin(L*h),np.cos(L*h),0,0],[0,0,1,L],[0,0,0,1]])
	else:
		matrix=np.array([[np.cos(L*KH),np.sin(L*KH)/KH,0,0],[-KH*np.sin(L*KH),np.cos(L*KH),0,0],[0,0,np.cos(L*(-K)**0.5),np.sin(L*(-K)**0.5)/(-K)**0.5],[0,0,-((-K)**0.5*np.sin(L*(-K)**0.5)),np.cos(L*(-K)**0.5)]])
	return matrix

def driff5(L):
	matrix= np.array([[1,L,0,0,0],[0,1,0,0,0],[0,0,1,L,0],[0,0,0,1,0],[0,0,0,0,1]])
	return matrix

def q5(L,K):
	if K>0:
		matrix= np.array([[math.cos(L*math.sqrt(K)),math.sin(L*math.sqrt(K))/math.sqrt(K),0,0,0],[-math.sqrt(K)*math.sin(L*math.sqrt(K)),math.cos(L*math.sqrt(K)),0,0,0],[0,0,math.cosh(L*math.sqrt(K)),math.sinh(L*math.sqrt(K))/math.sqrt(K),0],[0,0,math.sqrt(K)*math.sinh(L*math.sqrt(K)),math.cosh(L*math.sqrt(K)),0],[0,0,0,0,1]])
	elif K<0:
		K=-K
		matrix= np.array([[math.cosh(L*math.sqrt(K)),math.sinh(L*math.sqrt(K))/math.sqrt(K),0,0,0],[math.sqrt(K)*math.sinh(L*math.sqrt(K)),math.cosh(L*math.sqrt(K)),0,0,0],[0,0,math.cos(L*math.sqrt(K)),math.sin(L*math.sqrt(K))/math.sqrt(K),0],[0,0,-math.sqrt(K)*math.sin(L*math.sqrt(K)),math.cos(L*math.sqrt(K)),0],[0,0,0,0,1]])
	else:
		matrix= np.array([[1,L,0,0,0],[0,1,0,0,0],[0,0,1,L,0],[0,0,0,1,0],[0,0,0,0,1]])
	return matrix

def cfsd5(L,K,T0,L0):
	h=(np.pi*T0)/(180*L0)
	KH=(K+h**2)**0.5
	if K==0:
		matrix=np.array([[np.cos(L*h),np.sin(L*h)/h,0,0,(1-np.cos(L*h))/h],[-h*np.sin(L*h),np.cos(L*h),0,0,np.sin(L*h)],[0,0,1,L,0],[0,0,0,1,0],[0,0,0,0,1]])
	else:
		matrix=np.array([[np.cos(L*KH),np.sin(L*KH)/KH,0,0,h*(1-np.cos(L*KH))/KH**2],[-KH*np.sin(L*KH),np.cos(L*KH),0,0,h*np.sin(L*KH)/KH],[0,0,np.cos(L*(-K)**0.5),np.sin(L*(-K)**0.5)/(-K)**0.5,0],[0,0,-((-K)**0.5*np.sin(L*(-K)**0.5)),np.cos(L*(-K)**0.5),0],[0,0,0,0,1]])
	return matrix



############## Courant-Snyder Tranport #########################
def mcsft(M):
	matrix=np.array([[M[0,0]**2,-2*M[0,0]*M[0,1],M[0,1]**2,0,0,0],[-M[0,0]*M[1,0],M[0,0]*M[1,1]+M[0,1]*M[1,0],-M[0,1]*M[1,1],0,0,0],[M[1,0]**2,-2*M[1,0]*M[1,1],M[1,1]**2,0,0,0],[0,0,0,M[2,2]**2,-2*M[2,2]*M[2,3],M[2,3]**2],[0,0,0,-M[2,2]*M[3,2],M[2,2]*M[3,3]+M[2,3]*M[3,2],-M[2,3]*M[3,3]],[0,0,0,M[3,2]**2,-2*M[3,2]*M[3,3],M[3,3]**2]])
	return matrix

def linfunc(cell,csfunc,dispfunc,linpart):
	slen=0.0
	phadx=0.0
	phady=0.0
	chromx=0.0
	chromy=0.0
	radintegral=np.array([0.0,0.0,0.0,0.0,0.0,0.0])
	csft=[]
	dispt=[]
	st=[]
	phadxt=[]
	phadyt=[]
	csft.append(csfunc)
	dispt.append(dispfunc)
	st.append(slen)
	phadxt.append(phadx)
	phadyt.append(phady)
	for i in cell:
		if cell[i].length>=linpart and cell[i].type!='multipole':
			nsteps=int(cell[i].length*1.0001//linpart)
			for j in range(0,nsteps):
				bxi=csfunc[0]
				byi=csfunc[3]
				if cell[i].K!=0:
					chromx=chromx+linpart*bxi*cell[i].K
					chromy=chromy+linpart*byi*cell[i].K
				if cell[i].T0!=0:
					radintegral[0]=radintegral[0]+linpart*dispfunc[0]*(cell[i].T0*np.pi/(180.0*cell[i].length))
					radintegral[1]=radintegral[1]+linpart*(cell[i].T0*np.pi/(180.0*cell[i].length))**2
					radintegral[2]=radintegral[2]+linpart*abs(cell[i].T0*np.pi/(180.0*cell[i].length))**3
					radintegral[3]=radintegral[3]+linpart*(1+2*cell[i].K*(180.0*cell[i].length/(cell[i].T0*np.pi))**2)*dispfunc[0]*(cell[i].T0*np.pi/(180.0*cell[i].length))**3
					radintegral[4]=radintegral[4]+linpart*(bxi*dispfunc[1]**2+2.0*csfunc[1]*dispfunc[0]*dispfunc[1]+csfunc[2]*dispfunc[0]**2)*abs(cell[i].T0*np.pi/(180.0*cell[i].length))**3
					radintegral[5]=radintegral[5]+linpart*(dispfunc[0]*cell[i].K)**2
				mt=cell[i].MS
				csfunc=cell[i].MT.dot(csfunc)
				dispfunc=cell[i].M5T.dot(dispfunc)
				csft.append(csfunc)
				dispt.append(dispfunc)
				slen=slen+linpart
				st.append(slen)
				bxf=csfunc[0]
				byf=csfunc[3]
				phadx=phadx+np.arcsin(mt[0,1]/(bxi*bxf)**0.5)
				phady=phady+np.arcsin(mt[2,3]/(byi*byf)**0.5)
				phadxt.append(phadx)
				phadyt.append(phady)
			if cell[i].length>linpart*nsteps:
				bxi=csfunc[0]
				byi=csfunc[3]
				if cell[i].K!=0:
					chromx=chromx+(cell[i].length-nsteps*linpart)*bxi*cell[i].K
					chromy=chromy+(cell[i].length-nsteps*linpart)*byi*cell[i].K
				if cell[i].T0!=0:
					radintegral[0]=radintegral[0]+(cell[i].length-nsteps*linpart)*dispfunc[0]*(cell[i].T0*np.pi/(180*cell[i].length))
					radintegral[1]=radintegral[1]+(cell[i].length-nsteps*linpart)*(cell[i].T0*np.pi/(180*cell[i].length))**2
					radintegral[2]=radintegral[2]+(cell[i].length-nsteps*linpart)*abs(cell[i].T0*np.pi/(180*cell[i].length))**3
					radintegral[3]=radintegral[3]+(cell[i].length-nsteps*linpart)*(1+2*cell[i].K*(180*cell[i].length/(cell[i].T0*np.pi))**2)*dispfunc[0]*(cell[i].T0*np.pi/(180*cell[i].length))**3
					radintegral[4]=radintegral[4]+(cell[i].length-nsteps*linpart)*(bxi*dispfunc[1]**2+2*csfunc[1]*dispfunc[0]*dispfunc[1]+csfunc[2]*dispfunc[0]**2)*abs(cell[i].T0*np.pi/(180*cell[i].length))**3
					radintegral[5]=radintegral[5]+(cell[i].length-nsteps*linpart)*(dispfunc[0]*cell[i].K)**2
				mt=cell[i].MSF
				csfunc=cell[i].MF.dot(csfunc)
				dispfunc=cell[i].M5F.dot(dispfunc)
				csft.append(csfunc)
				dispt.append(dispfunc)
				slen=slen+cell[i].length-nsteps*linpart
				st.append(slen)
				bxf=csfunc[0]
				byf=csfunc[3]
				phadx=phadx+np.arcsin(mt[0,1]/(bxi*bxf)**0.5)
				phady=phady+np.arcsin(mt[2,3]/(byi*byf)**0.5)
				phadxt.append(phadx)
				phadyt.append(phady)
		if cell[i].length<linpart and cell[i].type!='multipole':
			bxi=csfunc[0]
			byi=csfunc[3]
			if cell[i].K!=0:
				chromx=chromx+cell[i].length*bxi*cell[i].K
				chromy=chromy+cell[i].length*byi*cell[i].K
			if cell[i].T0!=0:
				radintegral[0]=radintegral[0]+cell[i].length*dispfunc[0]*(cell[i].T0*np.pi/(180*cell[i].length))
				radintegral[1]=radintegral[1]+cell[i].length*(cell[i].T0*np.pi/(180*cell[i].length))**2
				radintegral[2]=radintegral[2]+cell[i].length*abs(cell[i].T0*np.pi/(180*cell[i].length))**3
				radintegral[3]=radintegral[3]+cell[i].length*(1+2*cell[i].K*(180*cell[i].length/(cell[i].T0*np.pi))**2)*dispfunc[0]*(cell[i].T0*np.pi/(180*cell[i].length))**3
				radintegral[4]=radintegral[4]+cell[i].length*(bxi*dispfunc[1]**2+2*csfunc[1]*dispfunc[0]*dispfunc[1]+csfunc[2]*dispfunc[0]**2)*abs(cell[i].T0*np.pi/(180*cell[i].length))**3
				radintegral[5]=radintegral[5]+cell[i].length*(dispfunc[0]*cell[i].K)**2
			mt=cell[i].MSF
			csfunc=cell[i].MF.dot(csfunc)
			dispfunc=cell[i].M5F.dot(dispfunc)
			csft.append(csfunc)
			dispt.append(dispfunc)
			slen=slen+cell[i].length
			st.append(slen)
			bxf=csfunc[0]
			byf=csfunc[3]
			phadx=phadx+np.arcsin(mt[0,1]/(bxi*bxf)**0.5)
			phady=phady+np.arcsin(mt[2,3]/(byi*byf)**0.5)
			phadxt.append(phadx)
			phadyt.append(phady)
		if cell[i].type=='multipole':
			bxi=csfunc[0]
			byi=csfunc[3]
			mt=cell[i].MS
			csfunc=cell[i].MT.dot(csfunc)
			dispfunc=cell[i].M5T.dot(dispfunc)
			csft.append(csfunc)
			dispt.append(dispfunc)
			slen=slen+cell[i].length
			st.append(slen)
			bxf=csfunc[0]
			byf=csfunc[3]
			phadx=phadx+np.arcsin(mt[0,1]/(bxi*bxf)**0.5)
			phady=phady+np.arcsin(mt[2,3]/(byi*byf)**0.5)
			phadxt.append(phadx)
			phadyt.append(phady)
	return phadx,phady,st,csft,dispt,radintegral,chromx,chromy


################ Correcion Cromatica ####################
def chrsext(schr1,schr2,chrx,chry,cell,csft,dispt,linpart,ntel,ncel):
    
	cintegral=np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
	idc=0
	for i in cell:
		if cell[i].type=='sextupole' or cell[i].K!=0:
			if cell[i].length>=linpart:
				nsteps=int(cell[i].length*1.0001//linpart)
				for j in range(0,nsteps):
					if cell[i].iloc[0]!=schr1 and cell[i].iloc[0]!=schr2:
						cintegral[0]=cintegral[0]+linpart*csft[idc][0]*(cell[i].K-2*dispt[idc][0]*cell[i].S)
						cintegral[1]=cintegral[1]+linpart*csft[idc][3]*(cell[i].K-2*dispt[idc][0]*cell[i].S)
					if cell[i].iloc[0]==schr1 or cell[i].iloc[0]==schr2:
						cintegral[2]=cintegral[2]+linpart*csft[idc][0]*(cell[i].K)
						cintegral[3]=cintegral[3]+linpart*csft[idc][3]*(cell[i].K)
					if cell[i].iloc[0]==schr1:
						cintegral[4]=cintegral[4]+2*linpart*csft[idc][0]*dispt[idc][0]
						cintegral[5]=cintegral[5]+2*linpart*csft[idc][3]*dispt[idc][0]
					if cell[i].iloc[0]==schr2:
						cintegral[6]=cintegral[6]+2*linpart*csft[idc][0]*dispt[idc][0]
						cintegral[7]=cintegral[7]+2*linpart*csft[idc][3]*dispt[idc][0]
					idc=idc+1
				if cell[i].length>linpart*nsteps:
					if cell[i].iloc[0]!=schr1 and cell[i].iloc[0]!=schr2:
						cintegral[0]=cintegral[0]+(cell[i].length-nsteps*linpart)*csft[idc][0]*(cell[i].K-2*dispt[idc][0]*cell[i].S)
						cintegral[1]=cintegral[1]+(cell[i].length-nsteps*linpart)*csft[idc][3]*(cell[i].K-2*dispt[idc][0]*cell[i].S)
					if cell[i].iloc[0]==schr1 or cell[i].iloc[0]==schr2:
						cintegral[2]=cintegral[2]+(cell[i].length-nsteps*linpart)*csft[idc][0]*(cell[i].K)
						cintegral[3]=cintegral[3]+(cell[i].length-nsteps*linpart)*csft[idc][3]*(cell[i].K)
					if cell[i].iloc[0]==schr1:
						cintegral[4]=cintegral[4]+2*(cell[i].length-nsteps*linpart)*csft[idc][0]*dispt[idc][0]
						cintegral[5]=cintegral[5]+2*(cell[i].length-nsteps*linpart)*csft[idc][3]*dispt[idc][0]
					if cell[i].iloc[0]==schr2:
						cintegral[6]=cintegral[6]+2*(cell[i].length-nsteps*linpart)*csft[idc][0]*dispt[idc][0]
						cintegral[7]=cintegral[7]+2*(cell[i].length-nsteps*linpart)*csft[idc][3]*dispt[idc][0]
					idc=idc+1
			if cell[i].length<linpart:
				if cell[i].iloc[0]!=schr1 and cell[i].iloc[0]!=schr2:
					cintegral[0]=cintegral[0]+cell[i].length*csft[idc][0]*(cell[i].K-2*dispt[idc][0]*cell[i].S)
					cintegral[1]=cintegral[1]+cell[i].length*csft[idc][3]*(cell[i].K-2*dispt[idc][0]*cell[i].S)
				if cell[i].iloc[0]==schr1 or cell[i].iloc[0]==schr2:
					cintegral[2]=cintegral[2]+cell[i].length*csft[idc][0]*(cell[i].K)
					cintegral[3]=cintegral[3]+cell[i].length*csft[idc][3]*(cell[i].K)
				if cell[i].iloc[0]==schr1:
					cintegral[4]=cintegral[4]+2*cell[i].length*csft[idc][0]*dispt[idc][0]
					cintegral[5]=cintegral[5]+2*cell[i].length*csft[idc][3]*dispt[idc][0]
				if cell[i].iloc[0]==schr2:
					cintegral[6]=cintegral[6]+2*cell[i].length*csft[idc][0]*dispt[idc][0]
					cintegral[7]=cintegral[7]+2*cell[i].length*csft[idc][3]*dispt[idc][0]
				idc=idc+1
		else:
			nsteps=int(cell[i].length*1.0001//linpart)
			if cell[i].length>=linpart:
				idc=idc+int(cell[i].length*1.0001//linpart)
				if cell[i].length>linpart*nsteps:
					idc=idc+1
			if cell[i].length<linpart:
				idc=idc+1
	cintegral[2]=cintegral[2]+cintegral[0]+chrx*(4*np.pi/int(ntel//ncel))
	cintegral[3]=cintegral[3]+cintegral[1]+chry*(-4*np.pi/int(ntel//ncel))
	schr1c=(cintegral[2]*cintegral[7]-cintegral[6]*cintegral[3])/(cintegral[4]*cintegral[7]-cintegral[6]*cintegral[5])
	schr2c=(cintegral[4]*cintegral[3]-cintegral[2]*cintegral[5])/(cintegral[4]*cintegral[7]-cintegral[6]*cintegral[5])
	return schr1c,schr2c

#### main function linear part ###
def build_linear_cell_and_functions(ring, elements, analysis, energy, linpart,
                                    whichprocess="nonlinear", plotlinear=0):

    if whichprocess == "nonlinear":
        global csftP, disptP, bxP, stP

    cell = pd.DataFrame()
    nn = 0

    for i in ring:
        if i[0][0] == 'ring':
            ntel = len(i[1])

        if i[0][0] == analysis:
            ncel = len(i[1])

            for j in i[1]:
                L = elements[j].iloc[1]
                rem = elements[j].length - int(elements[j].length * 1.0001 // linpart) * linpart

                if elements[j].type == 'drift':
                    cell[nn] = [j,'drift',L,0,0,0,0,0,0,
                                driff(L), driff5(L),
                                mcsft(driff(linpart)), mcsft(driff(rem)),
                                driff5(linpart), driff5(rem),
                                driff(linpart), driff(rem)]

                elif elements[j].type == 'bending':
                    T0 = elements[j].iloc[2]
                    K  = elements[j].iloc[3]
                    T1 = elements[j].iloc[4]
                    T2 = elements[j].iloc[5]
                    rem = L - int(L * 1.0001 // linpart) * linpart

                    cell[nn] = [j,'bending',L,T0,K,T1,T2,0,0,
                                cfsd(L,K,T0,L).real,
                                cfsd5(L,K,T0,L).real,
                                mcsft(cfsd(linpart,K,T0,L).real),
                                mcsft(cfsd(rem,K,T0,L).real),
                                cfsd5(linpart,K,T0,L).real,
                                cfsd5(rem,K,T0,L).real,
                                cfsd(linpart,K,T0,L).real,
                                cfsd(rem,K,T0,L).real]

                elif elements[j].type == 'quadrupole':
                    K = elements[j].iloc[3]

                    cell[nn] = [j,'quadrupole',L,0,K,0,0,0,0,
                                q(L,K), q5(L,K),
                                mcsft(q(linpart,K)), mcsft(q(rem,K)),
                                q5(linpart,K), q5(rem,K),
                                q(linpart,K), q(rem,K)]

                elif elements[j].type == 'sextupole':
                    S = elements[j].iloc[6]

                    cell[nn] = [j,'sextupole',L,0,0,0,0,S,0,
                                driff(L), driff5(L),
                                mcsft(driff(linpart)), mcsft(driff(rem)),
                                driff5(linpart), driff5(rem),
                                driff(linpart), driff(rem)]

                elif elements[j].type == 'multipole':
                    O = elements[j].iloc[7]

                    cell[nn] = [j,'multipole',L,0,0,0,0,0,O,
                                np.identity(4), np.identity(5),
                                mcsft(np.identity(4)), mcsft(np.identity(4)),
                                np.identity(5), np.identity(5),
                                np.identity(4), np.identity(4)]

                elif elements[j].type == 'octupole':
                    O = elements[j].iloc[7]

                    cell[nn] = [j,'multipole',L,0,0,0,0,0,O,
                                np.identity(4), np.identity(5),
                                mcsft(np.identity(4)), mcsft(np.identity(4)),
                                np.identity(5), np.identity(5),
                                np.identity(4), np.identity(4)]

                nn += 1

    cell.index = ['name','type','length','T0','K','T1','T2','S','O',
                  'M','M5','MT','MF','M5T','M5F','MS','MSF']

    m_cell = np.identity(4)
    m_cell5 = np.identity(5)

    for i in cell:
        m_cell = np.matmul(cell[i].M, m_cell)
        m_cell5 = np.matmul(cell[i].M5, m_cell5)

    sinmux = np.sign(m_cell[0,1]) * (
        -m_cell[0,1]*m_cell[1,0] - (m_cell[0,0]-m_cell[1,1])**2/4
    )**0.5

    sinmuy = np.sign(m_cell[2,3]) * (
        -m_cell[2,3]*m_cell[3,2] - (m_cell[2,2]-m_cell[3,3])**2/4
    )**0.5

    ax = (m_cell[0,0] - m_cell[1,1]) / (2*sinmux)
    ay = (m_cell[2,2] - m_cell[3,3]) / (2*sinmuy)

    bx = m_cell[0,1] / sinmux
    by = m_cell[2,3] / sinmuy

    gx = (1 + ax**2) / bx
    gy = (1 + ay**2) / by

    global beta_x_user, alpha_x_user, gamma_x_user
    beta_x_user = bx
    alpha_x_user = ax
    gamma_x_user = gx

    disp = (m_cell5[0,1]*m_cell5[1,4] + m_cell5[0,4]*(1-m_cell5[1,1])) / (
        2 - m_cell5[0,0] - m_cell5[1,1]
    )

    dispd = (m_cell5[1,0]*m_cell5[0,4] + m_cell5[1,4]*(1-m_cell5[0,0])) / (
        2 - m_cell5[0,0] - m_cell5[1,1]
    )

    csfunc = np.array([bx, ax, gx, by, ay, gy])
    dispfunc = np.array([disp, dispd, 0.0, 0.0, 1.0])

    phadx, phady, stP, csftP, disptP, radintegral, chromx, chromy = linfunc(
        cell, csfunc, dispfunc, linpart
    )

    bxt = extractelem(csftP, 0)
    byt = extractelem(csftP, 3)
    dist = [i * 100 for i in extractelem(disptP, 0)]

    Jx = 1 - radintegral[3] / radintegral[1]
    Js = 2 + radintegral[3] / radintegral[1]

    natemit = 3.8319e-13 * (1000*energy/0.5109989)**2 * radintegral[4] / (
        radintegral[1] - radintegral[3]
    )

    try:
        columns = os.get_terminal_size().columns
    except OSError:
        columns = 100

    print('=' * columns)
    print('Beta functions at s=0'.center(columns))
    print('| Ax =', ax, '| Ay =', ay, '| Bx =', bx, '| By =', by,
          '| Gx =', gx, '| Gy =', gy, '|')
    print('=' * columns)
    print('Ring'.center(columns))
    print('energy --- tunes --- chromaticities --- emittance --- circumference'.center(columns))
    print('| Energy =', energy,
          '|| Nux =', int(ntel//ncel)*phadx/(2*np.pi),
          '| Nuy =', int(ntel//ncel)*phady/(2*np.pi),
          '|| Chromx =', (-1/(4*np.pi))*int(ntel//ncel)*chromx,
          '| Chromy =', (1/(4*np.pi))*int(ntel//ncel)*chromy,
          '|| Emitx =', natemit,
          '|| circ =', int(ntel//ncel)*stP[len(stP)-1])
    print('=' * columns)

    if plotlinear == 1:
        plt.plot(stP, bxt, label='$\\beta_x$', linewidth=3.0, color='blue')
        plt.plot(stP, byt, label='$\\beta_y$', linewidth=3.0, color='red')
        plt.plot(stP, dist, label='$100*D$', linewidth=3.0, color='green')
        plt.xlabel('$s$ (m)', fontsize=20)
        plt.ylabel('Linear Functions (m)', fontsize=20)
        plt.legend()
        plt.show()

    return cell, csfunc, dispfunc, stP, csftP, disptP, bx, ax, gx, by, ay, gy, disp, dispd

#################################################################
#NEW FUNCTIONS CODE
##################################################################

######## Conversion: Polynomial <-> Vector + Hamiltonean Dict #######
        #Optional if dictionary is provided.
		
def poly_to_vector(poly, vars, vec_to_idx):
    if not isinstance(poly, sp.Poly):
        poly = sp.Poly(poly, *vars)

    vec = np.zeros(len(vec_to_idx), dtype=object)

    for monom, coeff in poly.terms():
        if monom in vec_to_idx:
            vec[vec_to_idx[monom]] = coeff

    return vec

def build_monomial_basis(vars, idx_to_vec):
    n = len(vars)
    basis = []

    for idx in range(len(idx_to_vec)):
        v = idx_to_vec[idx]
        term = 1
        for i in range(n):
            if v[i] != 0:
                term *= vars[i]**v[i]
        basis.append(term)

    return basis

def vector_to_poly(vec, monomial_basis):
    expr = 0

    for coeff, basis in zip(vec, monomial_basis):
        if coeff != 0:
            expr += coeff * basis

    return expr

def hamiltonean_dict(vars, H_pol, vec_to_idx):
  H_vec=poly_to_vector(H_pol, vars, vec_to_idx)
  H_dict = {}
  for i in range(len(H_vec)):
    if H_vec[i] != 0:
     H_dict[i] = H_vec[i]
  return H_dict

######### Principal index map functions #########
def indexmap(v):
    a = len(v)
    s = sum(v)
    index = 0

    for t in range(s):
        index += comb(t + a - 1, a - 1)

    remaining_sum = s

    for i in range(a - 1):
        for val in range(v[i] + 1, remaining_sum + 1):
            index += comb(remaining_sum - val + (a - i - 2), a - i - 2)
        remaining_sum -= v[i]

    return index

def vectormap(idx, a):

    s = 0
    while True:
        count = comb(s + a - 1, a - 1)
        if idx < count:
            break
        idx -= count
        s += 1

    v = []
    remaining_sum = s
    for i in range(a - 1):
        for val in range(remaining_sum, -1, -1):
            count = comb(remaining_sum - val + (a - i - 2), a - i - 2)
            if count <= idx:
                idx -= count
            else:
                v.append(val)
                remaining_sum -= val
                break

    v.append(remaining_sum)
    return v

################# General Non-linear Objects ############################

def load(m, d, H_pol, a_box, n=2):

    ndim = 2 * n

    idx_to_vec = {}
    vec_to_idx = {}

    layer_size = 0
    for s in range(m + 1):
        layer_size += comb(s + ndim - 1, ndim - 1)

    for x_delta in range(d + 1):
        offset = x_delta * layer_size

        for idx in range(layer_size):
            v = [x_delta] + vectormap(idx, ndim)
            v_tuple = tuple(v)

            global_idx = offset + idx
            idx_to_vec[global_idx] = v
            vec_to_idx[v_tuple] = global_idx

    size = len(idx_to_vec)
    dim = len(idx_to_vec[0])
	
  
    bracket_pairs = {}
    keys = list(idx_to_vec.keys())

    for i in keys:
        fi = idx_to_vec[i]

        for j in keys:
            fj = idx_to_vec[j]
            base = [fi[t] + fj[t] for t in range(dim)]

            for v in range(n):
                a_idx = 1 + v
                b_idx = 1 + n + v

                if base[a_idx] == 0 or base[b_idx] == 0:
                    continue

                cand = base.copy()
                cand[a_idx] -= 1
                cand[b_idx] -= 1

                k = vec_to_idx.get(tuple(cand))
                if k is not None:
                    bracket_pairs.setdefault((k, i), []).append((j, v))
					
    H_dict=hamiltonean_dict(vars, H_pol, vec_to_idx)

    G = np.zeros((size, size), dtype=float)

    for i in range(size):
        fi = idx_to_vec[i]
        for j in range(i, size):
            fj = idx_to_vec[j]

            val = 1.0
            for l in range(dim):
                s_ij = fi[l] + fj[l]
                if s_ij % 2 == 1:
                    val = 0.0
                    break
                val *= np.sqrt((2 * fi[l] + 1) * (2 * fj[l] + 1)) / (s_ij + 1)

            G[i, j] = val
            G[j, i] = val

    if np.isscalar(a_box):
        a_box = np.full(dim, float(a_box))
    else:
        a_box = np.array(a_box, dtype=float)

    C = np.zeros(size, dtype=float)

    for i in range(size):
        fi = idx_to_vec[i]
        val = 1.0
        for l in range(dim):
            val *= np.sqrt((2 * fi[l] + 1) / (2 * (a_box[l] ** (2 * fi[l] + 1))))
        C[i] = val
    C=C / C[0]

    HatC = [sps.lil_matrix((size, size), dtype=float) for _ in range(size)]

    for (k, i), pairs in bracket_pairs.items():
        fi = idx_to_vec[i]
        Ci_over_Ck = C[i] / C[k]

        for j, v in pairs:
            fj = idx_to_vec[j]

            a_idx = 1 + v
            b_idx = 1 + n + v

            sympl = fi[a_idx] * fj[b_idx] - fi[b_idx] * fj[a_idx]
            if sympl == 0:
                continue

            HatC[k][i, j] += sympl * Ci_over_Ck * C[j]

    HatC = [M.tocsr() for M in HatC]


    order = sorted(H_dict.keys())
    H_vec = sp.Matrix([sp.sympify(H_dict[j]) for j in order])

    HatC_T = [Mk.transpose().tocsr() for Mk in HatC]

    M_basis = []

    for j in order:
        Mj = sps.lil_matrix((size, size), dtype=float)

        for k in range(size):
            row = HatC_T[k].getrow(j)

            if row.nnz == 0:
                continue

            for i, val in zip(row.indices, row.data):
                Mj[k, i] = val / C[j]

        M_basis.append(Mj.tocsr())
    return idx_to_vec, vec_to_idx, bracket_pairs, G, C, HatC, H_vec, M_basis, order, H_dict


######### Construction of Non-LInear matrices ###########

def assemble_M(h_vec, M_basis):
	if len(M_basis) == 0:
		raise ValueError('M_basis is empty.')
	shape = M_basis[0].shape
	M = sps.csr_matrix(shape, dtype=float)
	for coeff, Mb in zip(h_vec, M_basis):
		if coeff != 0:
			M = M + float(coeff) * Mb
	return M.toarray()


def nulm2e(elem, H_vec, M_basis, tol=None, check_uper_right=False):
    if tol is None:
        tol = nulm2e_tol

    """
      [ Tqq   Tqn ]     # The block Tqn is zero.
      [ Tnq   Tnn ]

    """
    size = M_basis[0].shape[0]

    if elem.length == 0:
        tmatrix = np.eye(size)

    else:
        L = float(elem.length)

        b1n = float(elem.T0) * np.pi / (180.0 * L)
        b2n = float(elem.K)
        b3n = float(elem.S)
        b4n = float(elem.O)
        b5n = 0.0

        h_vec = np.array(H_vec_func(b1n, b2n, b3n, b4n, b5n), dtype=float).reshape(-1)

        M = assemble_M(h_vec, M_basis)
        tmatrix = scipy.linalg.expm(L * M)

    Mqq = tmatrix[:15, :15]
    Mqn = tmatrix[:15, 15:]
    Mnq = tmatrix[15:, :15]
    Mnn = tmatrix[15:, 15:]

    if check_uper_right and not np.all(np.abs(Mqn) < tol):
        print("Warning: el bloque superior derecho no es cero dentro de la tolerancia.")

    return tmatrix, Mqq, Mnn, Mqn, Mnq

def Non_linear_Transfer(cell, H_vec, M_basis):
	size = M_basis[0].shape[0]
	nl_cell = np.eye(size)
	for i in cell:
		tmatrix, _, _, _, _ = nulm2e(cell[i], H_vec, M_basis)
		nl_cell = np.dot(tmatrix, nl_cell)
	return nl_cell, nl_cell[15:,15:], nl_cell[15:,:15]

############ Objective Functions ################


def invariant(tnn, tnq, Sx, Sy, tol=None):
    if tol is None:
        tol = invariant_lstsq_tol

    Ux = tnq @ Sx
    Uy = tnq @ Sy

    D = np.eye(nonquad_size, dtype=float) - tnn

    # G = L L^T
    Lg = np.linalg.cholesky(Gnn)

    # ||U - D h||_G^2 = ||L^T (U - D h)||_2^2
    Aw = Lg.T @ D
    bx = Lg.T @ Ux
    by = Lg.T @ Uy

    hx, *_ = np.linalg.lstsq(Aw, bx, rcond=tol)
    hy, *_ = np.linalg.lstsq(Aw, by, rcond=tol)

    residual_x = Ux - D @ hx
    residual_y = Uy - D @ hy

    nhx = np.sqrt(max(hx @ Gnn @ hx, 0.0))
    nhy = np.sqrt(max(hy @ Gnn @ hy, 0.0))

    err_x = np.sqrt(max(residual_x.T @ Gnn @ residual_x, 0.0))
    err_y = np.sqrt(max(residual_y.T @ Gnn @ residual_y, 0.0))

    norm_Ux_G = np.sqrt(max(Ux.T @ Gnn @ Ux, 0.0))
    norm_Uy_G = np.sqrt(max(Uy.T @ Gnn @ Uy, 0.0))

    rel_err_x = err_x / max(norm_Ux_G, norm_floor_tol)
    rel_err_y = err_y / max(norm_Uy_G, norm_floor_tol)

    # Full invariant vectors:
    # Ix = quadratic Courant-Snyder part + nonlinear correction hx
    # Iy = quadratic Courant-Snyder part + nonlinear correction hy
    Ix = np.concatenate((Sx, hx))
    Iy = np.concatenate((Sy, hy))

    # Compute bracket_vec[k] = sum_{i,j} Ix[i] HatC[k][i,j] Iy[j]
    bracket_vec = np.zeros(len(Ix), dtype=float)

    for k in range(len(Ix)):
        if HatC[k].nnz != 0:
            bracket_vec[k] = Ix @ (HatC[k] @ Iy)

    bracket_sq = bracket_vec.T @ G @ bracket_vec

    if bracket_sq < 0.0 and abs(bracket_sq) < invariant_bracket_zero_tol:
        bracket_sq = 0.0

    n_bracketed = np.sqrt(max(bracket_sq, 0.0))

    return nhx, nhy, err_x, err_y, norm_Ux_G, norm_Uy_G, rel_err_x, rel_err_y, n_bracketed

def set_varied_parameters(v, change_quadrupoles=False):
    for j in range(len(varelem)):
        for i in cell:
            if cell[i].iloc[0] == varelem[j] and cell[i].type == 'sextupole':
                cell[i].S = v[j]

            if cell[i].iloc[0] == varelem[j] and cell[i].type == 'multipole':
                cell[i].O = v[j] * 1e8

            if change_quadrupoles and cell[i].iloc[0] == varelem[j] and cell[i].type == 'quadrupole':
                cell[i].K = v[j]
                cell[i].M = q(cell[i].length, cell[i].K)
                cell[i].M5 = q5(cell[i].length, cell[i].K)
                cell[i].MT = mcsft(q(linpart, cell[i].K))
                rem = cell[i].length - int(cell[i].length * 1.0001 // linpart) * linpart
                cell[i].MF = mcsft(q(rem, cell[i].K))
                cell[i].M5T = q5(linpart, cell[i].K)
                cell[i].M5F = q5(rem, cell[i].K)
                cell[i].MS = q(linpart, cell[i].K)
                cell[i].MSF = q(rem, cell[i].K)


def linear_data_for_objfunc():
    m_cell = np.identity(4)
    m_cell5 = np.identity(5)

    for i in cell:
        m_cell = np.matmul(cell[i].M, m_cell)
        m_cell5 = np.matmul(cell[i].M5, m_cell5)

    trace_x = m_cell[0, 0] + m_cell[1, 1]
    trace_y = m_cell[2, 2] + m_cell[3, 3]
    trace_x_excess = max(abs(trace_x) - max_linear_trace, 0.0)
    trace_y_excess = max(abs(trace_y) - max_linear_trace, 0.0)

    if trace_x_excess > 0.0 or trace_y_excess > 0.0:
        fobj = (
            trace_x_excess**2 * linear_trace_penalty_scale
            + trace_y_excess**2 * linear_trace_penalty_scale
        )
        return None, fobj

    sinmux = np.sign(m_cell[0, 1]) * (
        -m_cell[0, 1] * m_cell[1, 0]
        - (m_cell[0, 0] - m_cell[1, 1])**2 / 4
    )**0.5

    sinmuy = np.sign(m_cell[2, 3]) * (
        -m_cell[2, 3] * m_cell[3, 2]
        - (m_cell[2, 2] - m_cell[3, 3])**2 / 4
    )**0.5

    ax = (m_cell[0, 0] - m_cell[1, 1]) / (2 * sinmux)
    ay = (m_cell[2, 2] - m_cell[3, 3]) / (2 * sinmuy)

    bx = m_cell[0, 1] / sinmux
    by = m_cell[2, 3] / sinmuy

    gx = (1 + ax**2) / bx
    gy = (1 + ay**2) / by

    disp = (
        m_cell5[0, 1] * m_cell5[1, 4]
        + m_cell5[0, 4] * (1 - m_cell5[1, 1])
    ) / (2 - m_cell5[0, 0] - m_cell5[1, 1])

    dispd = (
        m_cell5[1, 0] * m_cell5[0, 4]
        + m_cell5[1, 4] * (1 - m_cell5[0, 0])
    ) / (2 - m_cell5[0, 0] - m_cell5[1, 1])

    csfunc = np.array([bx, ax, gx, by, ay, gy])
    dispfunc = np.array([disp, dispd, 0.0, 0.0, 1.0])

    phadx, phady, st, csft, dispt, radintegral, chromx, chromy = linfunc(
        cell, csfunc, dispfunc, linpart
    )

    byt = extractelem(csft, 3)
    dist = extractelem(dispt, 0)

    natemit = (
        3.8319e-13
        * (1000 * energy / 0.5109989)**2
        * radintegral[4]
        / (radintegral[1] - radintegral[3])
    )

    if (
        natemit > max_natemit
        or min(dist) < min_dispersion
        or dist[0] > max_initial_dispersion
        or max(byt) > max_beta_y
        or ax > max_alpha_x
    ):
        fobj = (
            natemit**2 * 1e41
            + abs(min(dist))**2 * 1e27
            + 2 * abs(ax)**2 * 1e23
            + abs(dist[0])**2 * 1e25
            + 0.3 * max(byt) * 1e20
        )
        return None, fobj

    return (bx, ax, gx, by, ay, gy, csft, dispt), None


def objfunc_given_linear_data(v, bx, ax, gx, by, ay, gy, csft, dispt,
                              print_table=False, save_csv=False,
                              return_details=False):
    try:
        term_size = os.get_terminal_size()
    except OSError:
        class TermSize:
            columns = 100
        term_size = TermSize()

    schr1c, schr2c = chrsext(
        schr1, schr2, chrx, chry, cell, csft, dispt, linpart, ntel, ncel
    )

    for i in cell:
        if cell[i].iloc[0] == schr1:
            cell[i].S = schr1c
        if cell[i].iloc[0] == schr2:
            cell[i].S = schr2c

    Sx = np.zeros(quad_size, dtype=float)
    Sx[idx_x2] = gx / C[idx_x2]
    Sx[idx_xpx] = 2.0 * ax / C[idx_xpx]
    Sx[idx_px2] = bx / C[idx_px2]

    Sy = np.zeros(quad_size, dtype=float)
    Sy[idx_y2] = gy / C[idx_y2]
    Sy[idx_ypy] = 2.0 * ay / C[idx_ypy]
    Sy[idx_py2] = by / C[idx_py2]

    Normx = np.sqrt(Sx.T @ Gqq @ Sx)
    Normy = np.sqrt(Sy.T @ Gqq @ Sy)

    _, tnn, tnq = Non_linear_Transfer(cell, H_vec, M_basis)

    nhx, nhy, err_x, err_y, norm_Ux_G, norm_Uy_G, rel_err_x, rel_err_y, n_bracketed = invariant(
        tnn, tnq, Sx, Sy
    )

    fobj = (
        (nhx + err_x * 10**2) / max(Normx, norm_floor_tol)
        + (nhy + err_y * 10**2) / max(Normy, norm_floor_tol)
    )

    details = {
        "schr1_name": schr1,
        "schr2_name": schr2,
        "schr1c": float(schr1c),
        "schr2c": float(schr2c),
        "nhx": float(nhx),
        "nhy": float(nhy),
        "err_x": float(err_x),
        "err_y": float(err_y),
        "norm_Ux_G": float(norm_Ux_G),
        "norm_Uy_G": float(norm_Uy_G),
        "rel_err_x": float(rel_err_x),
        "rel_err_y": float(rel_err_y),
        "Normx": float(Normx),
        "Normy": float(Normy),
        "n_bracketed": float(n_bracketed),
        "fobj": float(fobj),
        "bx": float(bx),
        "ax": float(ax),
        "gx": float(gx),
        "by": float(by),
        "ay": float(ay),
        "gy": float(gy),
    }

    if print_table:
        print('=' * term_size.columns)
        print('Chromatic correction'.center(term_size.columns))
        print('=' * term_size.columns)
        print(f"{schr1:<20} {schr1c:>25.16e}")
        print(f"{schr2:<20} {schr2c:>25.16e}")

        print('=' * term_size.columns)
        print('Nonlinear quantities'.center(term_size.columns))
        print('=' * term_size.columns)
        print(f"{'Quantity':<20} {'Value':>25}")
        print("-" * 45)
        print(f"{'nhx':<20} {nhx:>25.16e}")
        print(f"{'nhy':<20} {nhy:>25.16e}")
        print(f"{'err_x':<20} {err_x:>25.16e}")
        print(f"{'err_y':<20} {err_y:>25.16e}")
        print(f"{'norm_Ux_G':<20} {norm_Ux_G:>25.16e}")
        print(f"{'norm_Uy_G':<20} {norm_Uy_G:>25.16e}")
        print(f"{'rel_err_x':<20} {rel_err_x:>25.16e}")
        print(f"{'rel_err_y':<20} {rel_err_y:>25.16e}")
        print(f"{'Normx':<20} {Normx:>25.16e}")
        print(f"{'Normy':<20} {Normy:>25.16e}")
        print(f"{'n_bracketed':<20} {n_bracketed:>25.16e}")
        print("-" * 45)
        print(f"{'fobj':<20} {fobj:>25.16e}")
        print('=' * term_size.columns)

    if save_csv:
        save_f_obj_dinamica_csv(
            "f_obj_dinamica.csv",
            v,
            fobj,
            details,
            generation=None,
            solution_fitness=None,
            solution_idx=None,
            stage="direct_objfunc_call"
        )

    if return_details:
        return fobj, details

    return fobj


def objfunc(v, print_table=False, save_csv=False, return_details=False):
    v = np.array(v, dtype=float)

    if whichprocess == 'linear':
        set_varied_parameters(v, change_quadrupoles=True)

        lin_data, penalty = linear_data_for_objfunc()
        if penalty is not None:
            if return_details:
                details = {"penalty": True, "fobj": float(penalty)}
                return penalty, details
            return penalty

    elif whichprocess == 'nonlinear':
        set_varied_parameters(v, change_quadrupoles=False)

        lin_data = (bxP, axP, gxP, byP, ayP, gyP, csftP, disptP)

    else:
        raise ValueError("whichprocess must be 'linear' or 'nonlinear'.")

    return objfunc_given_linear_data(
        v,
        *lin_data,
        print_table=print_table,
        save_csv=save_csv,
        return_details=return_details
    )


########## Optimization ##################


import sys
import shutil


def _term_columns():
    return shutil.get_terminal_size((100, 20)).columns


def fitness_func(ga_instance, solution, solution_idx):
    fobj = objfunc(np.array(solution, dtype=float), print_table=False)

    fitness = 1.0 / (1.0 + fobj)

    gen_now = ga_instance.generations_completed + 1
    gen_total = ga_instance.num_generations
    pop_size = ga_instance.sol_per_pop
    calc_now = (solution_idx % pop_size) + 1
    left_now = pop_size - calc_now

    print(
        f"\rGeneration {gen_now}/{gen_total} | "
        f"calculation {calc_now}/{pop_size} | "
        f"left {left_now} | "
        f"f_obj = {fobj:.6e}",
        end="",
        flush=True
    )

    return fitness


def print_full_solution(label, solution):
    solution = np.array(solution, dtype=float)

    print(label, flush=True)
    print(
        np.array2string(
            solution,
            precision=17,
            separator=', ',
            suppress_small=False,
            max_line_width=100000
        ),
        flush=True
    )


def _write_dataframe_row(csv_name, row):
    df = pd.DataFrame([row])
    file_exists = os.path.exists(csv_name)

    df.to_csv(
        csv_name,
        mode="a",
        index=False,
        header=not file_exists,
        float_format="%.17e"
    )


def _solution_to_named_values(solution, prefix=""):
    solution = np.array(solution, dtype=float)
    row = {}

    for name, value in zip(varelem, solution):
        row[prefix + name] = float(value)

    return row


def save_resultados_csv(csv_name, solution, fobj, generation=None,
                        solution_fitness=None, solution_idx=None, stage="generation"):
    """
    Guarda el vector v0/best-solution al final de cada generación.
    Este archivo es el historial de soluciones, no el historial de métricas internas.
    """
    row = {
        "stage": stage,
        "generation": generation if generation is not None else "",
    }

    row.update(_solution_to_named_values(solution))
    row["fobj"] = float(fobj)

    if solution_fitness is not None:
        row["fitness"] = float(solution_fitness)

    if solution_idx is not None:
        row["solution_idx"] = int(solution_idx)

    _write_dataframe_row(csv_name, row)
    print(f"Saved resultados CSV in: {os.path.abspath(csv_name)}", flush=True)


# Alias for compatibility with your previous calls.
def save_solution_csv(csv_name, solution, fobj, solution_fitness=None, solution_idx=None, stage="generation"):
    save_resultados_csv(
        csv_name,
        solution,
        fobj,
        generation=None,
        solution_fitness=solution_fitness,
        solution_idx=solution_idx,
        stage=stage
    )


def save_f_obj_dinamica_csv(csv_name, solution, fobj, details, generation=None,
                            solution_fitness=None, solution_idx=None, stage="generation"):
    """
    Guarda las cantidades internas que aparecen en la tabla de objfunc/invariant
    para la mejor solución de cada generación.
    """
    row = {
        "stage": stage,
        "generation": generation if generation is not None else "",
    }

    if solution_fitness is not None:
        row["fitness"] = float(solution_fitness)

    if solution_idx is not None:
        row["solution_idx"] = int(solution_idx)

    row.update(_solution_to_named_values(solution, prefix="v_"))

    if details is not None:
        for key, value in details.items():
            if isinstance(value, (int, float, np.integer, np.floating)):
                row[key] = float(value)
            else:
                row[key] = value

    row["fobj"] = float(fobj)

    _write_dataframe_row(csv_name, row)
    print(f"Saved f_obj dinamica CSV in: {os.path.abspath(csv_name)}", flush=True)


def _element_type_from_cell(element_name):
    for i in cell:
        if cell[i].iloc[0] == element_name:
            return cell[i].type
    return None


def _opa_parameter_name(element_name):
    return "K" + str(element_name).upper()


def _opa_physical_value(element_name, optimizer_value):
    """
    Convierte el valor interno del optimizador al parametro que aparece
    en el archivo OPA original.

    En tu OPA local:
      - Sextupoles:  K = KD*3      => KD = valor_interno / 3.
      - Octupoles:   K = KO*1e-8   => KO = valor_interno * 1e8.
      - Quadrupoles: K = KQ        => KQ = valor_interno.
    """
    elem_type = _element_type_from_cell(element_name)
    value = float(optimizer_value)

    if elem_type == "multipole":
        return value * 1e8

    if elem_type == "sextupole":
        return value / 3.0

    return value


def write_opa_block_txt(txt_name, solution, details=None):
    """
    Escribe un bloque listo para pegar en el archivo .opa original.
    Se sobreescribe cada vez para que siempre tengas el último mejor bloque.
    """
    solution = np.array(solution, dtype=float)

    if len(solution) != len(varelem):
        raise ValueError(
            "The solution vector length does not match varelem. "
            f"len(solution)={len(solution)}, len(varelem)={len(varelem)}"
        )

    lines = []
    lines.append("{-----Parametros optimizados-----}")
    lines.append(f"{{-----Proceso: {whichprocess}-----}}")
    lines.append("{-----Variables optimizadas-----}")

    for name, value in zip(varelem, solution):
        pname = _opa_parameter_name(name)
        pvalue = _opa_physical_value(name, value)
        lines.append(f"{pname} ={pvalue:.16e};")

    if details is not None:
        optimized_names = set(varelem)
        has_schr1 = "schr1c" in details and schr1 not in optimized_names
        has_schr2 = "schr2c" in details and schr2 not in optimized_names

        if has_schr1 or has_schr2:
            lines.append("{-----Correccion cromatica-----}")

        if has_schr1:
            lines.append(f"{_opa_parameter_name(schr1)} ={_opa_physical_value(schr1, details['schr1c']):.16e};")
        if has_schr2:
            lines.append(f"{_opa_parameter_name(schr2)} ={_opa_physical_value(schr2, details['schr2c']):.16e};")

    lines.append("{-----Fin parametros optimizados-----}")

    with open(txt_name, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    print(f"Saved OPA block in: {os.path.abspath(txt_name)}", flush=True)


def _reset_output_files(file_names):
    for file_name in file_names:
        if file_name is not None and os.path.exists(file_name):
            os.remove(file_name)


def optimize(v0, num_generations=5, sol_per_pop=20, num_parents_mating=3,
             save_csv=True,
             csv_name="resultados.csv",
             fobj_csv_name="f_obj_dinamica.csv",
             opa_txt_name="bloque_opa_final.txt",
             reset_saved_files=True):

    v0 = np.array(v0, dtype=float)
    nvars = len(v0)

    if save_csv and reset_saved_files:
        _reset_output_files([csv_name, fobj_csv_name, opa_txt_name])

    v1 = [v0.copy()]

    for i in range(sol_per_pop - 1):
        v2 = np.zeros(nvars, dtype=float)

        if i == 0:
            for j in range(nvars):
                v2[j] = (1.0 + 0.0*np.random.normal()) * v0[j]

        elif i < int(0.6 * sol_per_pop):
            for j in range(nvars):
                v2[j] = (1.0 + 0.15*np.random.normal()) * v0[j]

        else:
            for j in range(nvars):
                v2[j] = (1.0 + 0.08*np.random.normal()) * v0[j]

        v1.append(v2)

    term_columns = _term_columns()

    print('=' * term_columns, flush=True)
    print('INITIAL VALUES'.center(term_columns), flush=True)
    print('=' * term_columns, flush=True)

    fobj0, details0 = objfunc(v0, print_table=True, return_details=True)

    print('=' * term_columns, flush=True)
    print('Initial solution summary'.center(term_columns), flush=True)
    print('=' * term_columns, flush=True)
    print(f"{'Initial f_obj':<30} {fobj0:>25.16e}", flush=True)
    print(f"{'Number of variables':<30} {nvars:>25}", flush=True)
    print('=' * term_columns, flush=True)

    print_full_solution("Initial solution vector:", v0)

    if save_csv:
        save_resultados_csv(
            csv_name,
            v0,
            fobj0,
            generation=0,
            solution_fitness=None,
            solution_idx=None,
            stage="initial"
        )
        save_f_obj_dinamica_csv(
            fobj_csv_name,
            v0,
            fobj0,
            details0,
            generation=0,
            solution_fitness=None,
            solution_idx=None,
            stage="initial"
        )
        write_opa_block_txt(opa_txt_name, v0, details0)

    def on_generation_callback(ga_instance):
        print("\n", flush=True)

        generation = ga_instance.generations_completed

        solution, solution_fitness, solution_idx = ga_instance.best_solution()
        solution = np.array(solution, dtype=float)

        fobj, details = objfunc(solution, print_table=True, return_details=True)

        term_columns = _term_columns()

        print('=' * term_columns, flush=True)
        print(('END OF GENERATION ' + str(generation)).center(term_columns), flush=True)
        print('=' * term_columns, flush=True)

        print(f"{'Generation number':<30} {generation:>25}", flush=True)
        print(f"{'Best f_obj':<30} {fobj:>25.16e}", flush=True)
        print(f"{'Best fitness':<30} {solution_fitness:>25.16e}", flush=True)
        print(f"{'Best solution index':<30} {solution_idx:>25}", flush=True)

        print('=' * term_columns, flush=True)
        print('Genetic algorithm summary'.center(term_columns), flush=True)
        print('=' * term_columns, flush=True)
        print(f"{'Population size':<30} {ga_instance.sol_per_pop:>25}", flush=True)
        print(f"{'Number of genes':<30} {ga_instance.num_genes:>25}", flush=True)
        print(f"{'Parents mating':<30} {ga_instance.num_parents_mating:>25}", flush=True)
        print(f"{'Selection type':<30} {ga_instance.parent_selection_type:>25}", flush=True)
        print(f"{'Crossover type':<30} {ga_instance.crossover_type:>25}", flush=True)
        print(f"{'Mutation type':<30} {ga_instance.mutation_type:>25}", flush=True)
        print(f"{'Mutation percent genes':<30} {ga_instance.mutation_percent_genes:>25}", flush=True)

        print('=' * term_columns, flush=True)
        print('Best solution vector'.center(term_columns), flush=True)
        print('=' * term_columns, flush=True)

        print_full_solution("Best solution vector:", solution)

        if save_csv:
            save_resultados_csv(
                csv_name,
                solution,
                fobj,
                generation=generation,
                solution_fitness=solution_fitness,
                solution_idx=solution_idx,
                stage="generation"
            )
            save_f_obj_dinamica_csv(
                fobj_csv_name,
                solution,
                fobj,
                details,
                generation=generation,
                solution_fitness=solution_fitness,
                solution_idx=solution_idx,
                stage="generation"
            )
            write_opa_block_txt(opa_txt_name, solution, details)

        print('=' * term_columns, flush=True)

    ga_instance = pygad.GA(
        num_generations=num_generations,
        num_parents_mating=num_parents_mating,
        fitness_func=fitness_func,
        sol_per_pop=sol_per_pop,
        num_genes=nvars,
        initial_population=np.array(v1, dtype=float),
        parent_selection_type="sss",
        crossover_type="single_point",
        mutation_type="random",
        mutation_percent_genes=20,
        keep_elitism=2,
        on_generation=on_generation_callback
    )

    ga_instance.run()

    print("\n", flush=True)

    solution, solution_fitness, solution_idx = ga_instance.best_solution()
    solution = np.array(solution, dtype=float)

    fobj, details = objfunc(solution, print_table=True, return_details=True)

    if fobj > fobj0:
        print("WARNING: GA result is worse than initial v0. Returning v0.", flush=True)
        solution = v0.copy()
        fobj = fobj0
        details = details0
        solution_fitness = 1.0 / (1.0 + fobj0)
        solution_idx = -1

    if save_csv:
        save_resultados_csv(
            csv_name,
            solution,
            fobj,
            generation=ga_instance.generations_completed,
            solution_fitness=solution_fitness,
            solution_idx=solution_idx,
            stage="final"
        )
        save_f_obj_dinamica_csv(
            fobj_csv_name,
            solution,
            fobj,
            details,
            generation=ga_instance.generations_completed,
            solution_fitness=solution_fitness,
            solution_idx=solution_idx,
            stage="final"
        )
        write_opa_block_txt(opa_txt_name, solution, details)

    term_columns = _term_columns()

    print('=' * term_columns, flush=True)
    print('FINAL BEST SOLUTION'.center(term_columns), flush=True)
    print('=' * term_columns, flush=True)
    print(f"{'Best fitness':<30} {solution_fitness:>25.16e}", flush=True)
    print(f"{'Best objective f_obj':<30} {fobj:>25.16e}", flush=True)
    print(f"{'Best solution index':<30} {solution_idx:>25}", flush=True)
    print('=' * term_columns, flush=True)

    print_full_solution("Final solution vector:", solution)

    print('=' * term_columns, flush=True)

    return solution, fobj, ga_instance

########################################################
#EXECUTABLE FUNCTIONS
########################################################3

def main():
	global cell, linpart, ntel, ncel, term_size, energy, chrx, chry, schr1, schr2, whichprocess, varelem, varelemplus, cellfile, madxfile
	global stP, csftP, disptP, bxP, axP, gxP, byP, ayP, gyP
	global idx_to_vec, vec_to_idx, bracket_pairs, G, C, HatC, H_vec, M_basis, order, H_dict, Gnn, Gqq
	global H_vec_func
	global idx_x2, idx_xpx, idx_px2, idx_y2, idx_ypy, idx_py2, quad_size, nonquad_size
	global nulm2e_tol, invariant_lstsq_tol, invariant_bracket_zero_tol, norm_floor_tol
	global max_natemit, min_dispersion, max_initial_dispersion, max_beta_y, max_alpha_x
	global max_linear_trace, linear_trace_penalty_scale

	# Basis settings
	m = 8  # Maximum polynomial degree/order used in the basis.
	d = 1  # Minimum order offset used when building the basis.
	n = 2  # Fixed for now.
	a_box = np.array([3e-2, 5e-3, 1.5e-3, 0.8e-3, 1e-3], dtype=float)  # Normalization/scaling box for phase-space variables.

	# File settings
	opafile = "opa_X_04-Feb-2025-10h53m_renglon_1.txt"  # Input OPA lattice file.
	# opafile = "esrf.txt"
	cellfile = "esrf.madx"  # Output MAD-X cell file name.
	madxfile = "trackingone"  # Base name used for tracking/MAD-X output.

	# Process settings
	analysis = "cell"  # Ring/cell label to analyze from the OPA file.
	linpart = 0.01  # Longitudinal step size for linear optics integration.
	whichprocess = "nonlinear"  # Selects linear or nonlinear optimization variables.
	plotlinear = 0  # Set to 1 to show linear optics plots.

	# Chromaticity settings
	chrx = 0  # Target horizontal chromaticity.
	chry = 0  # Target vertical chromaticity.
	schr1 = "d2"  # First sextupole family used for chromatic correction.
	schr2 = "d3"  # Second sextupole family used for chromatic correction.

	# Tolerance and limit settings
	nulm2e_tol = 1e-12  # Threshold for treating nonlinear matrix terms as zero.
	invariant_lstsq_tol = 1e-12  # Least-squares tolerance for invariant solving.
	invariant_bracket_zero_tol = 1e-25  # Tolerance for zero Poisson-bracket residual.
	norm_floor_tol = 1e-300  # Lower bound to avoid division by zero in norms.

	max_natemit = 71e-12  # Maximum allowed natural emittance.
	min_dispersion = 0.0  # Minimum allowed dispersion.
	max_initial_dispersion = 1e-2  # Maximum allowed starting dispersion.
	max_beta_y = 22  # Maximum allowed vertical beta function.
	max_alpha_x = 1e-3  # Maximum allowed horizontal alpha magnitude.
	max_linear_trace = 1.9  # Maximum allowed linear trace before penalty.
	linear_trace_penalty_scale = 1e28  # Penalty scale for linear trace violations.

	# Element settings
	base_varelem = [
		"qf1", "qd2", "qf3",
		"qf1i", "qd2i", "qf3i",
		"qf1d", "qd2d", "qf3d",
		"d1", "d4", "d3d", "d4d",
		"o2", "o3", "o4", "o3d", "o4d",
	]

	linear_varelem = [
		"qf1", "qd2", "qf3",
		"qf1i", "qd2i", "qf3i",
		"qf1d", "qd2d", "qf3d",
		"d1", "d4", "d3d", "d4d",
		"o2", "o3", "o4", "o3d", "o4d",
	]

	nonlinear_varelem = [
		"d1", "d4", "d3d", "d4d",
		"o2", "o3", "o4", "o3d", "o4d",
	]

	extra_varelemplus = [
		"momentum",
		"xamplitude",
	]

	# Initial optimization vectors
	v0_linear = [
		3.742153457895532e+00,
		-7.933030565890353e+00,
		9.383357074512245e+00,
		-1.961256869785227e+00,
		-7.622311324145445e+00,
		9.502339759326459e+00,
		1.528008695772636e+00,
		-6.400798167795065e+00,
		1.049016046874180e+01,
		3.148346595634115e+03 * 3,
		-4.621730496569649e+02 * 3,
		1.366797632793281e+02 * 3,
		-1.713849194855147e+02 * 3,
		-2.526041097299107e+12 * 1e-8,
		4.591682297815766e+07 * 1e-8,
		-3.136093044673090e+11 * 1e-8,
		-7.400633479124564e+10 * 1e-8,
		6.886521769314347e+11 * 1e-8,
	]

	v0_nonlinear = [
		3.148346595634115e+03 * 3,
		-4.621730496569649e+02 * 3,
		1.366797632793281e+02 * 3,
		-1.713849194855147e+02 * 3,
		-2.526041097299107e+12 * 1e-8,
		4.591682297815766e+07 * 1e-8,
		-3.136093044673090e+11 * 1e-8,
		-7.400633479124564e+10 * 1e-8,
		6.886521769314347e+11 * 1e-8,
	]

	# Optimization settings
	num_generations = 200  # Number of GA generations to run.
	sol_per_pop = 20  # Number of candidate solutions per generation.
	num_parents_mating = 4  # Number of parents selected for mating.
	save_csv = True  # Save optimization history files when True.
	csv_name = "resultados.csv"  # CSV with best solution values.
	fobj_csv_name = "f_obj_dinamica.csv"  # CSV with objective-function details.
	opa_txt_name = "bloque_opa_final.txt"  # Text block with final OPA parameters.
	reset_saved_files = True  # Delete previous output files before starting.

	H_pol = (
		sp.Rational(1, 2) * (px**2 + py**2) * (1 - delta + delta**2)
		- b1 * x * delta
		+ sp.Rational(1, 2) * b1**2 * x**2
		+ sp.Rational(1, 2) * b2 * (x**2 - y**2)
		+ sp.Rational(1, 3) * b3 * (x**3 - 3 * x * y**2)
		+ sp.Rational(1, 4) * b4 * (x**4 - 6 * x**2 * y**2 + y**4)
	)

	varelem = base_varelem
	varelemplus = varelem.copy()
	varelemplus.append(schr1)
	varelemplus.append(schr2)
	varelemplus = ['k' + item for item in varelemplus]
	varelemplus.extend(extra_varelemplus)

	if whichprocess == "linear":
		varelem = linear_varelem
		v0 = v0_linear
	else:
		varelem = nonlinear_varelem
		v0 = v0_nonlinear

	# loading
	ring, elements, energy = oparing(opafile)
	idx_to_vec, vec_to_idx, bracket_pairs, G, C, HatC, H_vec, M_basis, order, H_dict = load(
		m,
		d,
		H_pol,
		a_box,
		n=n,
	)

	H_vec_func = sp.lambdify((b1, b2, b3, b4, b5), H_vec, "numpy")

	for i in ring:
		if i[0][0] == 'ring':
			ntel = len(i[1])
		if i[0][0] == analysis:
			ncel = len(i[1])

	# linear
	cell, csfunc, dispfunc, stP, csftP, disptP, bx, ax, gx, by, ay, gy, disp, dispd = build_linear_cell_and_functions(
		ring,
		elements,
		analysis,
		energy,
		linpart,
		whichprocess,
		plotlinear=plotlinear,
	)

	bxP = bx
	axP = ax
	gxP = gx
	byP = by
	ayP = ay
	gyP = gy

	# Si son globales aparentemente no tienes que recuperarlas.
	Gnn = G[15:, 15:]
	Gqq = G[:15, :15]

	idx_x2 = vec_to_idx[(0, 2, 0, 0, 0)]
	idx_xpx = vec_to_idx[(0, 1, 0, 1, 0)]
	idx_px2 = vec_to_idx[(0, 0, 0, 2, 0)]

	idx_y2 = vec_to_idx[(0, 0, 2, 0, 0)]
	idx_ypy = vec_to_idx[(0, 0, 1, 0, 1)]
	idx_py2 = vec_to_idx[(0, 0, 0, 0, 2)]

	quad_size = 15
	nonquad_size = len(idx_to_vec) - 15

	# fobj = objfunc(v0, print_table=True)
	# quad_non_quad(bx, ax, gx, by, ay, gy)
	# print("Fobj:")
	# print(fobj)

	# quad_non_quad(bx, ax, gx, by, ay, gy)

	solution, fobj, ga_instance = optimize(
		v0,
		num_generations=num_generations,
		sol_per_pop=sol_per_pop,
		num_parents_mating=num_parents_mating,
		save_csv=save_csv,
		csv_name=csv_name,
		fobj_csv_name=fobj_csv_name,
		opa_txt_name=opa_txt_name,
		reset_saved_files=reset_saved_files,
	)

	print("Final solution: ")
	print(solution)
	print("Final fobj:")
	print(fobj)
	# solution, fobj, ga_instance = optimize(v0, 5, 20)

	# print("Final solution: ")
	# print(solution)
	

########################################################

############## Testing Functions ###############3333
def quad_non_quad(bx, ax, gx, by, ay, gy):
    t, tnn, tnq = Non_linear_Transfer(cell, H_vec, M_basis)
    tqq = t[:quad_size, :quad_size]

    Sx = np.zeros(quad_size, dtype=float)
    Sx[idx_x2]  = gx / C[idx_x2]
    Sx[idx_xpx] = 2.0 * ax / C[idx_xpx]
    Sx[idx_px2] = bx / C[idx_px2]

    Sy = np.zeros(quad_size, dtype=float)
    Sy[idx_y2]  = gy / C[idx_y2]    #no normalize courant snyder
    Sy[idx_ypy] = 2.0 * ay / C[idx_ypy]
    Sy[idx_py2] = by / C[idx_py2]

    rx = Sx - tqq @ Sx
    ry = Sy - tqq @ Sy

    nhx, nhy, err_x, err_y, *_ = invariant(tnn, tnq, Sx, Sy)

    print("Sx - tqq @ Sx =")
    print(rx)

    print("Sy - tqq @ Sy =")
    print(ry)

    print("norm x =", np.linalg.norm(rx))
    print("norm y =", np.linalg.norm(ry))

    print("Traza:")
    print(np.trace(tqq[5:,5:]))
    print("=====Non_quad:=========")
    print(nhx,nhy,err_x,err_y)

#exec(open("LS3.py").read())

#Tqn is the one that is zero.

if __name__=="__main__":
	main()
	
