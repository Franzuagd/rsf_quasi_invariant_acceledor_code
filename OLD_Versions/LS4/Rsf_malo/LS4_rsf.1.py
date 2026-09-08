######################## IMPORTS #####################################

import subprocess
import re
import os

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
###VEry wrong
#### main function linear part ###


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
                a_idx = 1+v
                b_idx = v+3

                if base[a_idx] == 0 or base[b_idx] == 0:
                    continue

                cand = base.copy()
                cand[a_idx] -= 1
                cand[b_idx] -= 1

                k = vec_to_idx.get(tuple(cand))
                if k is not None:
                    bracket_pairs.setdefault((k,i), []).append((j, v))

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

    epsilon = np.zeros(size, dtype=float)

    for i in range(size):
        fi = idx_to_vec[i]
        val = 1.0
        for l in range(dim):
            val *= np.sqrt((2 * fi[l] + 1) / (a_box[l] ** (2 * fi[l])))
        epsilon[i] = val

    B = [sps.lil_matrix((size, size), dtype=float) for _ in range(size)]

    for (k,i), pairs in bracket_pairs.items():
        fi = idx_to_vec[i]

        for j, v in pairs:
            fj = idx_to_vec[j]

            a_idx = 1+v
            b_idx = 3 + v

            sympl = fi[a_idx] * fj[b_idx] - fi[b_idx] * fj[a_idx]
            if sympl == 0:
                continue

            B[k][i, j] += sympl * epsilon[i]*epsilon[j]/epsilon[k]

    B = [M.tocsr() for M in B]


    order = sorted(H_dict.keys())
    H_vec = sp.Matrix([sp.sympify(H_dict[j]) for j in order])
    H_vec_func = sp.lambdify([b1, b2, b3, b4, b5], H_vec, modules='numpy')

    M_basis = []
    for i in order:
        Mi = sps.lil_matrix((size, size), dtype=float)

        for k in range(size):
            row = B[k].getrow(i)

            if row.nnz == 0:
                continue

            for j, val in zip(row.indices, row.data):
                Mi[k, j] = val/epsilon[i]

        M_basis.append(Mi.toarray())
	
    return idx_to_vec, vec_to_idx, bracket_pairs, G, epsilon, B, H_vec_func, M_basis, order, H_dict


######### Construction of Non-LInear matrices ###########

def assemble_M(h_vec, M_basis):
    if len(M_basis) == 0:
        raise ValueError('M_basis is empty.')

    shape = M_basis[0].shape
    M = np.zeros(shape, dtype=float)

    for coeff, Mb in zip(h_vec, M_basis):
        if coeff != 0:
            M += float(coeff) * Mb

    return M


def nulm2e(elem, H_vec_func, M_basis, tol=1e-14, check_uper_right=False):
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
        print("Warning: el bloque superior derecho no es cero.")

    return tmatrix, Mqq, Mnn, Mqn, Mnq

def Non_linear_Transfer(cell, H_vec, M_basis):
	size = M_basis[0].shape[0]
	nl_cell = np.eye(size)
	for i in cell:
		tmatrix, _, _, _, _ = nulm2e(cell[i], H_vec, M_basis)
		nl_cell = np.dot(tmatrix, nl_cell)
	return nl_cell, nl_cell[15:,15:], nl_cell[15:,:15]

############ Objective Functions ################


def invariant(tnn, tnq, Sx, Sy, tol=1e-14):

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

    rel_err_x = err_x / max(norm_Ux_G, 1e-300)
    rel_err_y = err_y / max(norm_Uy_G, 1e-300)

    # Full invariant vectors:
    # Ix = quadratic Courant-Snyder part + nonlinear correction hx
    # Iy = quadratic Courant-Snyder part + nonlinear correction hy
    Ix = np.concatenate((Sx, hx))
    Iy = np.concatenate((Sy, hy))

    # Compute bracket_vec[k] = sum_{i,j} Ix[i] HatC[k][i,j] Iy[j]
    bracket_vec = np.zeros(len(Ix), dtype=float)

    for k in range(len(Ix)):
        if B[k].nnz != 0:
            bracket_vec[k] = Ix @ (B[k] @ Iy)

    bracket_sq = bracket_vec.T @ G @ bracket_vec

    if bracket_sq < 0.0 and abs(bracket_sq) < 1e-25:
        bracket_sq = 0.0

    n_bracketed = np.sqrt(max(bracket_sq, 0.0))

    return nhx, nhy, err_x, err_y, norm_Ux_G, norm_Uy_G, rel_err_x, rel_err_y, n_bracketed, Ix, Iy

def set_varied_parameters(v, change_quadrupoles=False):
    for j in range(len(varelem)):
        for i in cell:
            if cell[i].iloc[0] == varelem[j] and cell[i].type == 'sextupole':
                cell[i].S = v[j]

            if cell[i].iloc[0] == varelem[j] and cell[i].type == 'multipole':
                cell[i].O = v[j]

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


def objfunc(v,print_table=False,plot_mode=False, return_details=False, testing_mode=False):
	v = np.array(v, dtype=float)

	try:
		term_size = os.get_terminal_size()
	except OSError:
		class TermSize:
			columns = 100
		term_size = TermSize()

#######################varied parameters#######################
	if whichprocess == 'nonlinear':
		set_varied_parameters(v, change_quadrupoles=False)

		csft = csftP
		dispt = disptP
		bx, ax, gx, by, ay, gy = map(float, csftP[0])

	if whichprocess=='linear':
		for j in range(0,len(varelem)):
			for i in cell:
				if cell[i].iloc[0]==varelem[j] and cell[i].type=='sextupole':
					cell[i].S=v[j]
				if cell[i].iloc[0]==varelem[j] and cell[i].type=='multipole':
					cell[i].O=v[j]
				if cell[i].iloc[0]==varelem[j] and cell[i].type=='quadrupole':
					cell[i].K=v[j]
					cell[i].M=q(cell[i].length,cell[i].K)
					cell[i].M5=q5(cell[i].length,cell[i].K)
					cell[i].MT=mcsft(q(linpart,cell[i].K))
					cell[i].MF=mcsft(q(cell[i].length-(int(cell[i].length*1.0001//linpart)*linpart),cell[i].K))
					cell[i].M5T=q5(linpart,cell[i].K)
					cell[i].M5F=q5(cell[i].length-(int(cell[i].length*1.0001//linpart)*linpart),cell[i].K)
					cell[i].MS=q(linpart,cell[i].K)
					cell[i].MSF=q(cell[i].length-(int(cell[i].length*1.0001//linpart)*linpart),cell[i].K)

#######################linear functions#######################
		m_cell=np.identity(4)
		m_cell5=np.identity(5)
		mux=0.0
		muy=0.0

		for i in cell:
			m_cell=np.matmul(cell[i].M,m_cell)
			m_cell5=np.matmul(cell[i].M5,m_cell5)

		if abs(m_cell[0,0]+m_cell[1,1])>2 or abs(m_cell[2,2]+m_cell[3,3])>2:
			fobj=abs(m_cell[0,0]+m_cell[1,1])*1e28+abs(m_cell[2,2]+m_cell[3,3])*1e28
			if return_details:
				return fobj, {"penalty": True, "reason": "unstable_trace", "fobj": float(fobj)}
			return fobj

		sinmux=np.sign(m_cell[0,1])*(-m_cell[0,1]*m_cell[1,0]-(m_cell[0,0]-m_cell[1,1])**2/4)**0.5
		sinmuy=np.sign(m_cell[2,3])*(-m_cell[2,3]*m_cell[3,2]-(m_cell[2,2]-m_cell[3,3])**2/4)**0.5
		ax=(m_cell[0,0]-m_cell[1,1])/2/sinmux
		ay=(m_cell[2,2]-m_cell[3,3])/2/sinmuy
		bx=m_cell[0,1]/sinmux
		by=m_cell[2,3]/sinmuy
		gx=(1+ax**2)/bx
		gy=(1+ay**2)/by

		disp=(m_cell5[0,1]*m_cell5[1,4]+m_cell5[0,4]*(1-m_cell5[1,1]))/(2-m_cell5[0,0]-m_cell5[1,1])
		dispd=(m_cell5[1,0]*m_cell5[0,4]+m_cell5[1,4]*(1-m_cell5[0,0]))/(2-m_cell5[0,0]-m_cell5[1,1])
		csfunc=np.array([bx,ax,gx,by,ay,gy])
		dispfunc=np.array([disp,dispd,0.0,0.0,1.0])

#######################displaying linear functions#######################
		phadx,phady,st,csft,dispt,radintegral,chromx,chromy=linfunc(cell,csfunc,dispfunc,linpart)
		bxt=extractelem(csft,0)
		byt=extractelem(csft,3)
		dist=extractelem(dispt,0)
		dist = [i * 1 for i in dist]

		Jx=1-(radintegral[3]/radintegral[1])
		Js=2+(radintegral[3]/radintegral[1])
		natemit=3.8319e-13*(1000*energy/0.5109989)**2*radintegral[4]/(radintegral[1]-radintegral[3])

		if natemit>71e-12 or min(dist)<0 or dist[0]>1e-2 or max(byt)>22 or ax > 1e-3:
			fobj=(natemit**2*1e41+abs(min(dist))**2*1e27+2*abs(ax)**2*1e23+abs(dist[0])**2*1e25+0.3*max(byt)*1e20)
			if return_details:
				return fobj, {
					"penalty": True,
					"reason": "linear_constraints",
					"natemit": float(natemit),
					"min_dist": float(min(dist)),
					"dist0": float(dist[0]),
					"max_byt": float(max(byt)),
					"ax": float(ax),
					"fobj": float(fobj),
				}
			return fobj

		if print_table:
			print('=' * term_size.columns)
			print('Beta functions at s=0'.center(term_size.columns))
			print('| Ax =',ax,'| Ay =',ay,'| Bx =',bx,'| By =',by,'| Gx =',gx,'| Gy =',gy,'|')
			print('=' * term_size.columns)
			print('Ring'.center(term_size.columns))
			print('energy --- tunes --- chromaticities --- emittance --- circumference'.center(term_size.columns))
			print('| Energy =',energy,'|| Nux =',int(ntel//ncel)*phadx/(2*np.pi),'| Nuy =',int(ntel//ncel)*phady/(2*np.pi),'|| Chromx =',(-1/(4*np.pi))*int(ntel//ncel)*chromx,'| Chromy =',(1/(4*np.pi))*int(ntel//ncel)*chromy,'|| Emitx =',natemit,'|| circ =',int(ntel//ncel)*st[len(st)-1])
			print('=' * term_size.columns)

#######################chromatic sextupoles calculation#######################
	schr1c,schr2c=chrsext(schr1,schr2,chrx,chry,cell,csft,dispt,linpart,ntel,ncel)

	param2madx=[]
	param2madx=v.tolist()
	param2madx.append(schr1c)
	param2madx.append(schr2c)
	param2madx.append(0)
	param2madx.append(6e-3)

	for i in cell:
		if cell[i].iloc[0]==schr1:
			cell[i].S=schr1c
		if cell[i].iloc[0]==schr2:
			cell[i].S=schr2c

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

	_, tnn, tnq = Non_linear_Transfer(cell, H_vec_func, M_basis)

	nhx, nhy, err_x, err_y, norm_Ux_G, norm_Uy_G, rel_err_x, rel_err_y, n_bracketed, Ix, Iy = invariant(
		tnn, tnq, Sx, Sy, tol=1e-14
	)

	fobj = nhx+err_x
	if plot_mode:
		plot_invariant_section(Ix, Iy, idx_to_vec, C, a_box, plane="both", levels=25)
	if testing_mode:
		return Ix,Iy
		
        

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


	if return_details:
		return fobj, details
        
	

	return fobj

            

    
######### Test functions ###############
def eval_plane(I, idx_to_vec, eps, Q, P, plane="x", delta0=0.0):
    plane = plane.lower()
    q, p, frozen = (1, 3, (2, 4)) if plane == "x" else (2, 4, (1, 3))
    Z = np.zeros_like(Q, dtype=float)

    for k, c in enumerate(I):
        if c == 0:
            continue

        v = idx_to_vec[k]

        if v[frozen[0]] != 0 or v[frozen[1]] != 0:
            continue

        term = float(c) * float(eps[k])

        if v[0] != 0:
            term *= delta0 ** v[0]
        if term == 0:
            continue
        if v[q] != 0:
            term *= Q ** v[q]
        if v[p] != 0:
            term *= P ** v[p]

        Z += term

    return Z
from datetime import datetime

def plot_invariant_section(
    Ix, Iy, idx_to_vec, eps, a_box,
    plane="both",
    levels=25,
    n=350,
    rmin=0.06,
    rmax=0.95,
    delta0=0.0,
    folder="graficos_invariantes"
):

    plane = plane.lower()

    if plane == "both":
        outx = plot_invariant_section(Ix, Iy, idx_to_vec, eps, a_box,
                                     "x", levels, n, rmin, rmax, delta0, folder)
        outy = plot_invariant_section(Ix, Iy, idx_to_vec, eps, a_box,
                                     "y", levels, n, rmin, rmax, delta0, folder)
        return {"x": outx, "y": outy}

    # ---------- Always plot over OPA aperture ----------
    aopa = [0, 5e-3, 2e-3, 1e-3, 1e-3]

    if plane == "x":
        I = Ix
        qmax = float(aopa[1])
        pmax = float(aopa[3])
        qlab, plab = "x", "px"
        title = "Ix on y=py=0"

    elif plane == "y":
        I = Iy
        qmax = float(aopa[2])
        pmax = float(aopa[4])
        qlab, plab = "y", "py"
        title = "Iy on x=px=0"

    else:
        raise ValueError("plane must be 'x', 'y', or 'both'.")

    qvals = np.linspace(-qmax, qmax, n)
    pvals = np.linspace(-pmax, pmax, n)
    Q, P = np.meshgrid(qvals, pvals, indexing="xy")

    Z = eval_plane(I, idx_to_vec, eps, Q, P, plane, delta0)

    if isinstance(levels, int):

        levels *= 2      # Double the contour count

        theta = np.linspace(0.0, 2.0*np.pi, 300, endpoint=False)
        radii = np.linspace(rmin, rmax, levels)
        lev = []

        for r in radii:
            Qr = r*qmax*np.cos(theta)
            Pr = r*pmax*np.sin(theta)

            vals = eval_plane(I, idx_to_vec, eps, Qr, Pr, plane, delta0)
            vals = vals[np.isfinite(vals)]

            if len(vals):
                lev.append(float(np.median(vals)))

        lev = np.unique(np.round(np.sort(np.array(lev)), 14))

        zmin = np.nanmin(Z)
        zmax = np.nanmax(Z)

        lev = lev[(lev > zmin) & (lev < zmax)]

        if len(lev) < 2:
            lev = np.linspace(zmin, zmax, levels + 2)[1:-1]

    else:
        lev = np.array(levels, dtype=float)

    os.makedirs(folder, exist_ok=True)

    plt.figure(figsize=(8, 6))

    plt.contour(Q, P, Z, levels=lev)      # No contour labels

    plt.xlabel(qlab)
    plt.ylabel(plab)
    plt.title(title)

    plt.xlim(-qmax, qmax)
    plt.ylim(-pmax, pmax)

    plt.grid(True)
    plt.tight_layout()

    filename = f"invariant_section_{plane}.png"
    path = os.path.join(folder, datetime.now().strftime("%Y%m%d_%H%M%S")+filename)

    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()

    return Q, P, Z, lev, path


def testing(v, tol=1e-14):
    v = np.array(v, dtype=float)
    csft=csftP; dispt=disptP

    if whichprocess == 'nonlinear':

        set_varied_parameters(v, change_quadrupoles=False)

        csft = csftP
        dispt = disptP

        bx, ax, gx, by, ay, gy = map(float, csftP[0])


    if whichprocess == 'linear':
        for j in range(0, len(varelem)):
            for i in cell:
                if cell[i].iloc[0] == varelem[j] and cell[i].type == 'sextupole':
                    cell[i].S = v[j]
                if cell[i].iloc[0] == varelem[j] and cell[i].type == 'multipole':
                    cell[i].O = v[j]
                if cell[i].iloc[0] == varelem[j] and cell[i].type == 'quadrupole':
                    cell[i].K = v[j]
                    cell[i].M = q(cell[i].length, cell[i].K)
                    cell[i].M5 = q5(cell[i].length, cell[i].K)
                    cell[i].MT = mcsft(q(linpart, cell[i].K))
                    cell[i].MF = mcsft(
                        q(
                            cell[i].length - int(cell[i].length * 1.0001 // linpart) * linpart,
                            cell[i].K
                        )
                    )
                    cell[i].M5T = q5(linpart, cell[i].K)
                    cell[i].M5F = q5(
                        cell[i].length - int(cell[i].length * 1.0001 // linpart) * linpart,
                        cell[i].K
                    )
                    cell[i].MS = q(linpart, cell[i].K)
                    cell[i].MSF = q(
                        cell[i].length - int(cell[i].length * 1.0001 // linpart) * linpart,
                        cell[i].K
                    )

        m_cell = np.identity(4)
        m_cell5 = np.identity(5)

        for i in cell:
            m_cell = np.matmul(cell[i].M, m_cell)
            m_cell5 = np.matmul(cell[i].M5, m_cell5)

        if abs(m_cell[0, 0] + m_cell[1, 1]) > 2 or abs(m_cell[2, 2] + m_cell[3, 3]) > 2:
            fobj = abs(m_cell[0, 0] + m_cell[1, 1]) * 1e28
            fobj += abs(m_cell[2, 2] + m_cell[3, 3]) * 1e28
            return fobj

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
            cell,
            csfunc,
            dispfunc,
            linpart
        )

        bxt = extractelem(csft, 0)
        byt = extractelem(csft, 3)
        dist = extractelem(dispt, 0)
        dist = [i * 1 for i in dist]

        Jx = 1 - radintegral[3] / radintegral[1]
        Js = 2 + radintegral[3] / radintegral[1]

        natemit = (
            3.8319e-13
            * (1000 * energy / 0.5109989)**2
            * radintegral[4]
            / (radintegral[1] - radintegral[3])
        )

        if natemit > 71e-12 or min(dist) < 0 or dist[0] > 1e-2 or max(byt) > 22 or ax > 1e-3:
            print(natemit, min(dist), dist[0], max(byt), ax)

            fobj = (
                natemit**2 * 1e41
                + abs(min(dist))**2 * 1e27
                + 2 * abs(ax)**2 * 1e23
                + abs(dist[0])**2 * 1e25
                + 0.3 * max(byt) * 1e20
            )

            return fobj

        term_size = os.get_terminal_size()

        print('=' * term_size.columns)
        print('Beta functions at s=0'.center(term_size.columns))
        print('| Ax =', ax, '| Ay =', ay, '| Bx =', bx, '| By =', by, '| Gx =', gx, '| Gy =', gy, '|')
        print('=' * term_size.columns)
        print('Ring'.center(term_size.columns))
        print('energy --- tunes --- chromaticities --- emittance --- circumference'.center(term_size.columns))
        print(
            '| Energy =', energy,
            '|| Nux =', int(ntel // ncel) * phadx / (2 * np.pi),
            '| Nuy =', int(ntel // ncel) * phady / (2 * np.pi),
            '|| Chromx =', (-1 / (4 * np.pi)) * int(ntel // ncel) * chromx,
            '| Chromy =', (1 / (4 * np.pi)) * int(ntel // ncel) * chromy,
            '|| Emitx =', natemit,
            '|| circ =', int(ntel // ncel) * st[len(st) - 1]
        )
        print('=' * term_size.columns)

    schr1c, schr2c = chrsext(
        schr1,
        schr2,
        chrx,
        chry,
        cell,
        csft,
        dispt,
        linpart,
        ntel,
        ncel
    )

    param2madx = []
    param2madx = v.tolist()
    param2madx.append(schr1c)
    param2madx.append(schr2c)
    param2madx.append(0)
    param2madx.append(6e-3)

    term_size = os.get_terminal_size()

    for i in cell:
        if cell[i].iloc[0] == schr1:
            cell[i].S = schr1c
        if cell[i].iloc[0] == schr2:
            cell[i].S = schr2c

    print('Chromatic correction'.center(term_size.columns))
    print('| ', schr1, '  =', schr1c, '| ', schr2, ' =', schr2c, '|')
    print('=' * term_size.columns)

    t, tnn, tnq = Non_linear_Transfer(cell, H_vec_func, M_basis)
    tqq = t[:quad_size, :quad_size]

    Sx = np.zeros(quad_size, dtype=float)
    Sx[idx_x2] = gx / C[idx_x2]
    Sx[idx_xpx] = 2.0 * ax / C[idx_xpx]
    Sx[idx_px2] = bx / C[idx_px2]

    Sy = np.zeros(quad_size, dtype=float)
    Sy[idx_y2] = gy / C[idx_y2]
    Sy[idx_ypy] = 2.0 * ay / C[idx_ypy]
    Sy[idx_py2] = by / C[idx_py2]

    sx = tqq @ Sx
    sy = tqq @ Sy

    rx = Sx - sx
    ry = Sy - sy

    print("Refrence/Error for x")
    print(Sx @ Gqq @ Sx)
    print(rx @ Gqq @ rx)
    print(sx @ Gqq @ sx)

    print("Reference/Error for y")
    print(Sy @ Gqq @ Sy)
    print(ry @ Gqq @ ry)
    print(sy @ Gqq @ sy)

    print("Matrix")
    print(Sx)
    print(Sy)
    print("epsilon")

    Ux = tnq @ Sx
    Uy = tnq @ Sy

    D = np.eye(nonquad_size, dtype=float) - tnn

    Lg = np.linalg.cholesky(Gnn)

    Aw = Lg.T @ D
    bx_ls = Lg.T @ Ux
    by_ls = Lg.T @ Uy

    hx, *_ = np.linalg.lstsq(Aw, bx_ls, rcond=tol)
    hy, *_ = np.linalg.lstsq(Aw, by_ls, rcond=tol)

    residual_x = Ux - D @ hx
    residual_y = Uy - D @ hy

    nhx = np.sqrt(max(hx @ Gnn @ hx, 0.0))
    nhy = np.sqrt(max(hy @ Gnn @ hy, 0.0))

    err_x = np.sqrt(max(residual_x.T @ Gnn @ residual_x, 0.0))
    err_y = np.sqrt(max(residual_y.T @ Gnn @ residual_y, 0.0))

    norm_Ux_G = np.sqrt(max(Ux.T @ Gnn @ Ux, 0.0))
    norm_Uy_G = np.sqrt(max(Uy.T @ Gnn @ Uy, 0.0))

    rel_err_x = err_x / max(norm_Ux_G, 1e-300)
    rel_err_y = err_y / max(norm_Uy_G, 1e-300)

    Ix = np.concatenate((Sx, hx))
    Iy = np.concatenate((Sy, hy))

    bracket_vec = np.zeros(len(Ix), dtype=float)

    for k in range(len(Ix)):
        if B[k].nnz != 0:
            bracket_vec[k] = Ix @ (B[k] @ Iy)

    bracket_sq = bracket_vec.T @ G @ bracket_vec

    if bracket_sq < 0.0 and abs(bracket_sq) < 1e-25:
        bracket_sq = 0.0

    n_bracketed = np.sqrt(max(bracket_sq, 0.0))

    print("Quasinvariant X: (Norm/rel_error/h)")
    print(nhx)
    print(rel_err_x)

    print("Quasinvariant Y: (Norm/rel_error/h)")
    print(nhy)
    print(rel_err_y)

    print("Checking it out")

    ix = t @ Ix
    iy = t @ Iy

    rrx = Ix - ix
    rry = Iy - iy

    print("Refrence/Error for x")
    print(Ix @ G @ Ix)
    print(rrx @ G @ rrx)
    print(ix @ G @ ix)

    print("Reference/Error for y")
    print(Iy @ G @ Iy)
    print(rry @ G @ rry)
    print(iy @ G @ iy)

    plot_invariant_section(Ix, Iy, idx_to_vec, C, a_box, plane="both", levels=25)

########## Optimization ##################


import sys
import shutil


def _term_columns():
    return shutil.get_terminal_size((100, 20)).columns



def _element_type_from_cell(element_name):
    for i in cell:
        if cell[i].iloc[0] == element_name:
            return cell[i].type
    return None


def _opa_parameter_name(element_name):
    return "K" + str(element_name).upper()


def write_opa_block_txt(txt_name, solution, details=None):

    solution = np.asarray(solution, dtype=float)

    lines = []
    lines.append("{-----Parametros optimizados-----}")

    # write optimized variables correctly
    for name, value in zip(varelem, solution):
        pname = _opa_parameter_name(name)
        lines.append(f"{pname} = {value:.16e};")

    # write extra parameters from details safely
    if details is not None:

        if "schr1c" in details:
            lines.append(f"{_opa_parameter_name(schr1)} = {details['schr1c']:.16e};")

        if "schr2c" in details:
            lines.append(f"{_opa_parameter_name(schr2)} = {details['schr2c']:.16e};")

    lines.append("{-----Fin parametros optimizados-----}")

    with open(txt_name, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    print(f"Saved OPA block in: {os.path.abspath(txt_name)}", flush=True)


def _reset_output_files(file_names):
    for file_name in file_names:
        if file_name is not None and os.path.exists(file_name):
            os.remove(file_name)



import numpy as np
from scipy.optimize import minimize
import cma


def hybrid_optimize(
		v0,
		sigma=0.5,
		cma_iters=200,
		popsize=None,
		save_every=20,
		print_every=10,
		opa_txt_name="OPA_BLOQUE_HY",
		powell_iters=200,
		powell_maxfev=None  # <-- NUEVO: límite de evaluaciones
):
	v0 = np.asarray(v0, dtype=float)

	# ===== HISTORIAL =====
	history = []  # Lista para guardar [iteración, fase, best_f, eval_counter]
	history.append([0, "inicial", float('inf'), 0])

	fobj0, details0 = objfunc(v0, print_table=True, return_details=True)
	write_opa_block_txt(opa_txt_name, v0, details0)

	# CONTADOR DE EVALUACIONES
	eval_counter = 0

	def objective(x):
		nonlocal eval_counter
		eval_counter += 1

		# Print cada 10 evaluaciones para monitorear
		if eval_counter % 10 == 0:
			print(f"[fobj #{eval_counter}] evaluando...")

		fobj, _ = objfunc(
			x,
			print_table=False,
			return_details=True
		)
		return float(fobj)

	# Track best solution
	best_x = v0.copy()
	best_f = objective(v0)
	history.append([0, "inicial", best_f, eval_counter])

	# CMA-ES
	es = cma.CMAEvolutionStrategy(
		v0,
		sigma,
		{
			"popsize": popsize if popsize else 4 + int(3 * np.log(len(v0))),
			"verb_disp": 0,
		}
	)

	for i in range(cma_iters):

		solutions = es.ask()
		values = [objective(x) for x in solutions]
		es.tell(solutions, values)

		idx = np.argmin(values)
		if values[idx] < best_f:
			best_f = values[idx]
			best_x = solutions[idx].copy()
			# Guardar en historial CADA VEZ que mejora
			history.append([i + 1, "cma_mejora", best_f, eval_counter])

		if i % print_every == 0:
			print(f"[CMA-ES] iter {i}/{cma_iters} | best f = {best_f:.6e} | evals = {eval_counter}")

			history.append([i + 1, "cma_print", best_f, eval_counter])

		if i % save_every == 0:
			_, details = objfunc(best_x, return_details=True)
			write_opa_block_txt(opa_txt_name, best_x, details)

	x_cma = best_x.copy()
	print(f"Starting Powell (maxiter={powell_iters}, maxfev={powell_maxfev if powell_maxfev else 'ilimitado'})")
	history.append([cma_iters, "inicio_powell", best_f, eval_counter])

	# Configuración de Powell CON límite de evaluaciones
	powell_options = {
		"maxiter": powell_iters,
		"disp": False
	}


	if powell_maxfev is not None:
		powell_options["maxfev"] = powell_maxfev


	res = minimize(
		objective,
		x_cma,
		method="Powell",
		options=powell_options
	)

	# Guardar mejor de Powell durante la optimización
	# (Powell no da acceso fácil a su historial interno)
	# Guardamos el resultado final
	x_final = res.x
	f_final = res.fun

	history.append(["powell", "resultado_powell", f_final, eval_counter])

	if f_final > best_f:
		x_final = best_x
		f_final = best_f
		history.append(["powell", "powell_empeoro", f_final, eval_counter])
	else:
		history.append(["powell", "powell_mejoro", f_final, eval_counter])

	_, details = objfunc(x_final, return_details=True)
	write_opa_block_txt(opa_txt_name, x_final, details)

	# ===== IMPRIMIR HISTORIAL AL FINAL =====
	print("\n" + "=" * 60)
	print(" HISTORIAL COMPLETO DE OPTIMIZACIÓN")
	print("=" * 60)
	print(f"{'Iteración':<12} | {'Fase':<18} | {'Best f':<15} | {'Evals':<8}")
	print("-" * 60)

	for entry in history:
		iter_num, fase, f_val, evals = entry
		if isinstance(iter_num, str):
			# Para entradas con string (powell)
			print(f"{iter_num:<12} | {fase:<18} | {f_val:.6e}   | {evals:<8}")
		else:
			print(f"{iter_num:<12} | {fase:<18} | {f_val:.6e}   | {evals:<8}")

	print("-" * 60)
	print(f"{'FINAL':<12} | {'':<18} | {f_final:.6e}   | {eval_counter:<8}")
	print("=" * 60)

	# También imprimir el resumen original
	print("\n" + "=" * 50)
	print("HYBRID OPTIMIZATION COMPLETE")
	print("=" * 50)
	print(f"Best CMA-ES f   : {best_f:.6e}")
	print(f"Final Powell f  : {f_final:.6e}")
	print(f"Total fobj evals: {eval_counter}")
	print("=" * 50)

	return x_final, f_final, es, res, history  # <-- AHORA RETORNA HISTORIAL



########################################################
#EXECUTABLE FUNCTIONS
########################################################3

###Issue: There exist a strong dependency of I with the box, this shouldn't happened at all.
def main():
	global cell, linpart, ntel, ncel, term_size, energy, chrx, chry, schr1, schr2, whichprocess, varelem, varelemplus, cellfile, madxfile
	global idx_to_vec, vec_to_idx, bracket_pairs, G, epsilon, B, H_vec_func, M_basis, order, H_dict, C, a_box
    #user setting
	m = 8
	d = 0
	n = 2  #Fixed for now! 
	a_box = np.array([0.05e-2, 3.5e-3, 1e-3, 0.8e-3, 0.6e-3], dtype=float)
	#a_box = np.array([0.2e-2, 5e-3, 5e-3, 1e-3, 1e-3], dtype=float)
	#a_box = np.array([0.1e-2, 4e-3, 0.5e-3, 0.8e-3, 0.2e-3], dtype=float) #x y px py
	#a_box = np.array([1, 1, 1, 1, 1], dtype=float) #x y px py
	H_pol = (
    sp.Rational(1, 2)*(px**2 + py**2)*(1 - delta + delta**2)
    - b1*x*delta
    + sp.Rational(1, 2)*b1**2*x**2
    + sp.Rational(1, 2)*b2*(x**2 - y**2)
    + sp.Rational(1, 3)*b3*(x**3 - 3*x*y**2)
    + sp.Rational(1, 4)*b4*(x**4 - 6*x**2*y**2 + y**4)
	)
      
	#whichprocess='linear'
	opafile='esrf2026_2_1_1116_BAD.opa'
	#opafile='esrf.txt'
	cellfile='esrf.madx'
	madxfile='trackingone'
	chrx=0
	chry=0

	#Change 2 Actual name
	schr1='sf1'
	schr2='sd1'; varelem=['qf1','qd2','qd3','qf4','qd5','qf7']
	#varelem=['qf1','qd2','qf3','qf1i','qd2i','qf3i','qf1d','qd2d','qf3d','d1','d4','d3d','d4d','o2','o3','o4','o3d','o4d']
	varelemplus=[]
	varelemplus=varelem.copy()
	varelemplus.append(schr1)
	varelemplus.append(schr2)
	varelemplus = ['k' + item for item in varelemplus]
	varelemplus.append('momentum')
	varelemplus.append('xamplitude')
	analysis='cell'
	linpart=0.01
	whichprocess = 'nonlinear'
	v0=[3.633167514008421e+00,-4.258277621861492e+00,-2.690860661253351e+00,2.754505457254375e+00,-3.336431720192718e+00,-1.492176552197721e+00]

    #loading
	ring,elements,energy=oparing(opafile)
	idx_to_vec, vec_to_idx, bracket_pairs, G, epsilon, B, H_vec_func, M_basis, order, H_dict = load(m, d, H_pol, a_box, n=2)

    #linear
	if whichprocess=='nonlinear':
		global csftP,disptP,bxP,stP
#######################calculating lattice linear matrices#######################
	cell=pd.DataFrame()
	nn=0
	for i in ring:
		if i[0][0]=='ring':
			ntel=len(i[1])
		if i[0][0]==analysis:
			ncel=len(i[1])
			for j in i[1]:
				if elements[j].type=='drift':
					cell[nn]=[j,'drift',elements[j].iloc[1],0,0,0,0,0,0,driff(elements[j].iloc[1]),driff5(elements[j].iloc[1]),mcsft(driff(linpart)),mcsft(driff(elements[j].length-(int(elements[j].length*1.0001//linpart)*linpart))),driff5(linpart),driff5(elements[j].length-(int(elements[j].length*1.0001//linpart)*linpart)),driff(linpart),driff(elements[j].length-(int(elements[j].length*1.0001//linpart)*linpart))]
				elif elements[j].type=='bending':
					cell[nn]=[j,'bending',elements[j].iloc[1],elements[j].iloc[2],elements[j].iloc[3],elements[j].iloc[4],elements[j].iloc[5],0,0,cfsd(elements[j].iloc[1],elements[j].iloc[3],elements[j].iloc[2],elements[j].iloc[1]).real,cfsd5(elements[j].iloc[1],elements[j].iloc[3],elements[j].iloc[2],elements[j].iloc[1]).real,mcsft(cfsd(linpart,elements[j].iloc[3],elements[j].iloc[2],elements[j].iloc[1]).real),mcsft(cfsd(elements[j].iloc[1]-(int(elements[j].iloc[1]*1.0001//linpart)*linpart),elements[j].iloc[3],elements[j].iloc[2],elements[j].iloc[1]).real),cfsd5(linpart,elements[j].iloc[3],elements[j].iloc[2],elements[j].iloc[1]).real,cfsd5(elements[j].iloc[1]-(int(elements[j].iloc[1]*1.0001//linpart)*linpart),elements[j].iloc[3],elements[j].iloc[2],elements[j].iloc[1]).real,cfsd(linpart,elements[j].iloc[3],elements[j].iloc[2],elements[j].iloc[1]).real,cfsd(elements[j].iloc[1]-(int(elements[j].iloc[1]*1.0001//linpart)*linpart),elements[j].iloc[3],elements[j].iloc[2],elements[j].iloc[1]).real]
				elif elements[j].type=='quadrupole':
					cell[nn]=[j,'quadrupole',elements[j].iloc[1],0,elements[j].iloc[3],0,0,0,0,q(elements[j].iloc[1],elements[j].iloc[3]),q5(elements[j].iloc[1],elements[j].iloc[3]),mcsft(q(linpart,elements[j].iloc[3])),mcsft(q(elements[j].length-(int(elements[j].length*1.0001//linpart)*linpart),elements[j].iloc[3])),q5(linpart,elements[j].iloc[3]),q5(elements[j].length-(int(elements[j].length*1.0001//linpart)*linpart),elements[j].iloc[3]),q(linpart,elements[j].iloc[3]),q(elements[j].length-(int(elements[j].length*1.0001//linpart)*linpart),elements[j].iloc[3])]
				elif elements[j].type=='sextupole':
					cell[nn]=[j,'sextupole',elements[j].iloc[1],0,0,0,0,elements[j].iloc[6],0,driff(elements[j].iloc[1]),driff5(elements[j].iloc[1]),mcsft(driff(linpart)),mcsft(driff(elements[j].length-(int(elements[j].length*1.0001//linpart)*linpart))),driff5(linpart),driff5(elements[j].length-(int(elements[j].length*1.0001//linpart)*linpart)),driff(linpart),driff(elements[j].length-(int(elements[j].length*1.0001//linpart)*linpart))]
				elif elements[j].type=='multipole':
					cell[nn]=[j,'multipole',elements[j].iloc[1],0,0,0,0,0,elements[j].iloc[7],np.identity(4),np.identity(5),mcsft(np.identity(4)),mcsft(np.identity(4)),np.identity(5),np.identity(5),np.identity(4),np.identity(4)]
				elif elements[j].type=='octupole':
					cell[nn]=[j,'multipole',elements[j].iloc[1],0,0,0,0,0,elements[j].iloc[7],np.identity(4),np.identity(5),mcsft(np.identity(4)),mcsft(np.identity(4)),np.identity(5),np.identity(5),np.identity(4),np.identity(4)]
				nn=nn+1
	cell.index=['name','type','length','T0','K','T1','T2','S','O','M','M5','MT','MF','M5T','M5F','MS','MSF']
#######################linear functions#######################
	m_cell=np.identity(4)
	m_cell5=np.identity(5)
	mux=0.0
	muy=0.0
	for i in cell:
		m_cell=np.matmul(cell[i].M,m_cell)
		m_cell5=np.matmul(cell[i].M5,m_cell5)
	sinmux=np.sign(m_cell[0,1])*(-m_cell[0,1]*m_cell[1,0]-(m_cell[0,0]-m_cell[1,1])**2/4)**0.5
	sinmuy=np.sign(m_cell[2,3])*(-m_cell[2,3]*m_cell[3,2]-(m_cell[2,2]-m_cell[3,3])**2/4)**0.5
	ax=(m_cell[0,0]-m_cell[1,1])/2/sinmux
	ay=(m_cell[2,2]-m_cell[3,3])/2/sinmuy
	bx=m_cell[0,1]/sinmux
	by=m_cell[2,3]/sinmuy
	gx=(1+ax**2)/bx
	gy=(1+ay**2)/by
	disp=(m_cell5[0,1]*m_cell5[1,4]+m_cell5[0,4]*(1-m_cell5[1,1]))/(2-m_cell5[0,0]-m_cell5[1,1])
	dispd=(m_cell5[1,0]*m_cell5[0,4]+m_cell5[1,4]*(1-m_cell5[0,0]))/(2-m_cell5[0,0]-m_cell5[1,1])
	csfunc=np.array([bx,ax,gx,by,ay,gy])
	dispfunc=np.array([disp,dispd,0.0,0.0,1.0])
#	print(dispfunc)
#######################displaying linear functions#######################
	phadx,phady,stP,csftP,disptP,radintegral,chromx,chromy=linfunc(cell,csfunc,dispfunc,linpart)
	bxt=extractelem(csftP,0)
	byt=extractelem(csftP,3)
#	print(max(byt))
	dist=extractelem(disptP,0)
#dispt=extractelem(dispt,1)
	dist = [i * 100 for i in dist]
#dispt = [i * 100 for i in dispt]
	Jx=1-(radintegral[3]/radintegral[1])
	Js=2+(radintegral[3]/radintegral[1])
	natemit=3.8319e-13*(1000*energy/0.5109989)**2*radintegral[4]/(radintegral[1]-radintegral[3])
	term_size = os.get_terminal_size()
	print('=' * term_size.columns)
	print('Beta functions at s=0'.center(term_size.columns))
	print('| Ax =',ax,'| Ay =',ay,'| Bx =',bx,'| By =',by,'| Gx =',gx,'| Gy =',gy,'|')
	print('=' * term_size.columns)
	print('Ring'.center(term_size.columns))
	print('energy --- tunes --- chromaticities --- emittance --- circumference'.center(term_size.columns))
	print('| Energy =',energy,'|| Nux =',int(ntel//ncel)*phadx/(2*np.pi),'| Nuy =',int(ntel//ncel)*phady/(2*np.pi),'|| Chromx =',(-1/(4*np.pi))*int(ntel//ncel)*chromx,'| Chromy =',(1/(4*np.pi))*int(ntel//ncel)*chromy,'|| Emitx =',natemit,'|| circ =',int(ntel//ncel)*stP[len(stP)-1])
	print('=' * term_size.columns)
#######################plot linear functions#######################
	plotlinear=0
	if plotlinear==1:
		plt.plot(stP,bxt,label='$\\beta_x$',linewidth=3.0,color='blue')
		plt.plot(stP,byt,label='$\\beta_y$',linewidth=3.0,color='red')
		plt.plot(stP,dist,label='$100*D$',linewidth=3.0,color='green')
#		plt.plot(st,dispt,label='$100*dD/ds$',linewidth=3.0,color='yellow')
		plt.xlabel('$s$ (m)', fontsize=20)
		plt.ylabel('Linear Functions (m)', fontsize=20)
		plt.legend()
		plt.show()

	bxP = bx
	global idx_x2, idx_xpx, idx_px2, idx_y2, idx_ypy, idx_py2, quad_size, nonquad_size, Gnn, Gqq
    
	
	Gnn=G[15:,15:]
	Gqq=G[:15,:15]
	idx_x2  = vec_to_idx[(0, 2, 0, 0, 0)]
	idx_xpx = vec_to_idx[(0, 1, 0, 1, 0)]
	idx_px2 = vec_to_idx[(0, 0, 0, 2, 0)]

	idx_y2  = vec_to_idx[(0, 0, 2, 0, 0)]
	idx_ypy = vec_to_idx[(0, 0, 1, 0, 1)]
	idx_py2 = vec_to_idx[(0, 0, 0, 0, 2)]

	quad_size = 15
	nonquad_size = len(idx_to_vec)-15
	normalization_coef=np.sqrt(bx*gx/(epsilon[idx_px2]*epsilon[idx_x2])-ax**2/(epsilon[idx_xpx]**2)) #fixing quadratic part to be close to 1
	epsilon=epsilon*normalization_coef
	C=epsilon #Temporar para que funcione

	
    #################################
	#fobj = objfunc(v0, print_table=True) 
	Ix,Iy=objfunc(v0,plot_mode=True,testing_mode=True)
	coeff1 = np.asarray(Ix, dtype=float) * np.asarray(C, dtype=float)
	coeff2 = np.asarray(Iy, dtype=float) * np.asarray(C, dtype=float)
	#print("Ix")
	#print(coeff1)
	#print("Ix(Base)")
	#print(Ix)
	opti=False
	if opti:
		vfinal, f, es, res, history = hybrid_optimize(
            v0,
            sigma=0.5,
            cma_iters=2,#20,
            popsize=4,#2,
            save_every=1,#20,
            print_every=1,#0,
            powell_iters=2,#20,
            powell_maxfev=5#300
        )
		testing(vfinal)

	# Extraer solo las mejoras
	#mejoras = [h for h in history if h[1] in ['inicial', 'cma_mejora', 'powell_mejoro']]
	#print(mejoras)


	#print("Fobj:")
	#print(fobj)


	#print("Final solution: ")
	#print(solution)
	

########################################################


#exec(open("LS3.py").read())

#Tqn is the one that is zero.

if __name__=="__main__":
	main()



##### SERIUS Mistake some functions still use the convention of e^8 of octupoles.
###Observation, Sx usually is very small, we can make it larger and it should affect the optimization process,
#This by chosing a normalization constant A based on twis parameters.