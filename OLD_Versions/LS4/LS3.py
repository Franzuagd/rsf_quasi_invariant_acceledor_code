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

    C=C/C[0]


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


def nulm2e(elem, H_vec, M_basis, tol=1e-12, check_uper_right=True):
    """
      [ Mqq   Mqn ]     # The block Mqn is zero.
      [ Mnq   Mnn ]

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

        h_vec_sym = H_vec.subs({
            b1: b1n,
            b2: b2n,
            b3: b3n,
            b4: b4n,
            b5: b5n
        })

        h_vec = np.array([float(x) for x in h_vec_sym], dtype=float)

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


def invariant(tnn, tnq, Sx, Sy, tol=1e-30):

    Ux = tnq @ Sx
    Uy = tnq @ Sy

    D = np.eye(nonquad_size, dtype=float) - tnn
    A = D.T @ Gnn @ D

    bx = D.T @ Gnn @ Ux
    by = D.T @ Gnn @ Uy

    hx, *_ = np.linalg.lstsq(A, bx, rcond=None)
    hy, *_ = np.linalg.lstsq(A, by, rcond=None)

    residual_x = Ux + tnn @ hx - hx
    residual_y = Uy + tnn @ hy - hy

    nhx = np.sqrt(hx @ Gnn @ hx)
    nhy = np.sqrt(hy @ Gnn @ hy)
    err_x = np.sqrt(residual_x.T @ Gnn @ residual_x)
    err_y = np.sqrt(residual_y.T @ Gnn @ residual_y)

    norm_Ux_G = np.sqrt(Ux.T @ Gnn @ Ux)
    norm_Uy_G = np.sqrt(Uy.T @ Gnn @ Uy)

    return nhx, nhy, err_x, err_y

def objfunc(v, print_table=False, save_csv=False):
	#######################varied parameters#######################
	if whichprocess=='nonlinear':
		for j in range(0,len(varelem)):
			for i in cell:
				if cell[i].iloc[0]==varelem[j] and cell[i].type=='sextupole':
					cell[i].S=v[j]
				if cell[i].iloc[0]==varelem[j] and cell[i].type=='multipole':
					cell[i].O=v[j]*1e8
	if whichprocess=='linear':
		for j in range(0,len(varelem)):
			for i in cell:
				if cell[i].iloc[0]==varelem[j] and cell[i].type=='sextupole':
					cell[i].S=v[j]
				if cell[i].iloc[0]==varelem[j] and cell[i].type=='multipole':
					cell[i].O=v[j]*1e8
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
#		print(dispfunc)
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
			print(natemit,min(dist),dist[0],max(byt),ax)
			fobj=(natemit**2*1e41+abs(min(dist))**2*1e27+2*abs(ax)**2*1e23+abs(dist[0])**2*1e25+0.3*max(byt)*1e20)
			return fobj
		term_size = os.get_terminal_size()
		print('=' * term_size.columns)
		print('Beta functions at s=0'.center(term_size.columns))
		print('| Ax =',ax,'| Ay =',ay,'| Bx =',bx,'| By =',by,'| Gx =',gx,'| Gy =',gy,'|')
		print('=' * term_size.columns)
		print('Ring'.center(term_size.columns))
		print('energy --- tunes --- chromaticities --- emittance --- circumference'.center(term_size.columns))
		print('| Energy =',energy,'|| Nux =',int(ntel//ncel)*phadx/(2*np.pi),'| Nuy =',int(ntel//ncel)*phady/(2*np.pi),'|| Chromx =',(-1/(4*np.pi))*int(ntel//ncel)*chromx,'| Chromy =',(1/(4*np.pi))*int(ntel//ncel)*chromy,'|| Emitx =',natemit,'|| circ =',int(ntel//ncel)*st[len(st)-1])
		print('=' * term_size.columns)
#######################chromatic sextupoles calculation#######################
	schr1c,schr2c=chrsext(schr1,schr2,chrx,chry,cell,csft,dispt,linpart,ntel,ncel) #csftP,disptP
	param2madx=[]
	param2madx=list(v)
	param2madx.append(schr1c)
	param2madx.append(schr2c)
	param2madx.append(0)
	param2madx.append(6e-3)
#	print(param2madx)
#	param2madx[-1]=0 #
#	tracking=getTracking(varelemplus,param2madx,cellfile,madxfile) #
#	print(tracking)
	term_size = os.get_terminal_size()
	for i in cell:
		if cell[i].iloc[0]==schr1:
			cell[i].S=schr1c
		if cell[i].iloc[0]==schr2:
			cell[i].S=schr2c
	print('Chromatic correction'.center(term_size.columns))
	print('| ',schr1,'  =',schr1c,'| ',schr2,' =',schr2c,'|')
	print('=' * term_size.columns)
#######################plot linear functions#######################
	plotlinear=0
	bxt=extractelem(csftP,0)
	byt=extractelem(csftP,3)
	dist=extractelem(disptP,0)
	dist = [i * 100 for i in dist]
	if plotlinear==1:
		plt.plot(stP,bxt,label='$\\beta_x$',linewidth=3.0,color='blue')
		plt.plot(stP,byt,label='$\\beta_y$',linewidth=3.0,color='red')
		plt.plot(stP,dist,label='$100*D$',linewidth=3.0,color='green')
#		plt.plot(st,dispt,label='$100*dD/ds$',linewidth=3.0,color='yellow')
		plt.xlabel('$s$ (m)', fontsize=20)
		plt.ylabel('Linear Functions (m)', fontsize=20)
		plt.legend()
		plt.show()

 ################Printing ################

	if print_table:
		print('=' * term_size.columns)
		print('Constraint penalty'.center(term_size.columns))
		print('=' * term_size.columns)
		print(f"{'natemit':<20} {natemit:>25.16e}")
		print(f"{'min(dist)':<20} {min(dist):>25.16e}")
		print(f"{'dist[0]':<20} {dist[0]:>25.16e}")
		print(f"{'max(byt)':<20} {max(byt):>25.16e}")
		print(f"{'ax':<20} {ax:>25.16e}")
		print('=' * term_size.columns)

		print('=' * term_size.columns)
		print('Beta functions at s=0'.center(term_size.columns))
		print('=' * term_size.columns)
		print(f"{'Ax':<20} {ax:>25.16e}")
		print(f"{'Ay':<20} {ay:>25.16e}")
		print(f"{'Bx':<20} {bx:>25.16e}")
		print(f"{'By':<20} {by:>25.16e}")
		print(f"{'Gx':<20} {gx:>25.16e}")
		print(f"{'Gy':<20} {gy:>25.16e}")
		print('=' * term_size.columns)

		print('Ring'.center(term_size.columns))
		print('=' * term_size.columns)
		print(f"{'Energy':<20} {energy:>25.16e}")
		print(f"{'Nux':<20} {int(ntel//ncel)*phadx/(2*np.pi):>25.16e}")
		print(f"{'Nuy':<20} {int(ntel//ncel)*phady/(2*np.pi):>25.16e}")
		print(f"{'Chromx':<20} {(-1/(4*np.pi))*int(ntel//ncel)*chromx:>25.16e}")
		print(f"{'Chromy':<20} {(1/(4*np.pi))*int(ntel//ncel)*chromy:>25.16e}")
		print(f"{'Emitx':<20} {natemit:>25.16e}")
		print(f"{'Circumference':<20} {int(ntel//ncel)*st[len(st)-1]:>25.16e}")
		print('=' * term_size.columns)
	if print_table:
		print('Chromatic correction'.center(term_size.columns))
		print('=' * term_size.columns)
		print(f"{schr1:<20} {schr1c:>25.16e}")
		print(f"{schr2:<20} {schr2c:>25.16e}")
		print('=' * term_size.columns)

	plotlinear=0

	bxt=extractelem(csftP,0)
	byt=extractelem(csftP,3)
	dist=extractelem(disptP,0)
	dist = [i * 100 for i in dist]

	if plotlinear==1:
		plt.plot(stP,bxt,label='$\\beta_x$',linewidth=3.0,color='blue')
		plt.plot(stP,byt,label='$\\beta_y$',linewidth=3.0,color='red')
		plt.plot(stP,dist,label='$100*D$',linewidth=3.0,color='green')
		plt.xlabel('$s$ (m)', fontsize=20)
		plt.ylabel('Linear Functions (m)', fontsize=20)
		plt.legend()
		plt.show()

	Sx = np.zeros(quad_size, dtype=float)
	Sx[idx_x2]  = gx / C[idx_x2]
	Sx[idx_xpx] = 2.0 * ax / C[idx_xpx]
	Sx[idx_px2] = bx / C[idx_px2]

	Sy = np.zeros(quad_size, dtype=float)
	Sy[idx_y2]  = gy / C[idx_y2]    # We need to normalize it
	Sy[idx_ypy] = 2.0 * ay / C[idx_ypy]
	Sy[idx_py2] = by / C[idx_py2]


	_,tnn, tnq=Non_linear_Transfer(cell, H_vec, M_basis) 
	nhx, nhy, err_x, err_y= invariant(tnn, tnq, Sx, Sy, tol=1e-12)

	fobj=nhx+nhy+err_x+err_y

	fobj = nhx + nhy + err_x + err_y

	if print_table:
		print('=' * term_size.columns)
		print('Fluctuation index'.center(term_size.columns))
		print('=' * term_size.columns)

		print(f"{'Quantity':<20} {'Value':>25}")
		print("-" * 45)
		print(f"{'nhx':<20} {nhx:>25.16e}")
		print(f"{'nhy':<20} {nhy:>25.16e}")
		print(f"{'err_x':<20} {err_x:>25.16e}")
		print(f"{'err_y':<20} {err_y:>25.16e}")
		print("-" * 45)
		print(f"{'fobj':<20} {fobj:>25.16e}")
		print('=' * term_size.columns)
	
	if save_csv:
		solution = np.array(v)
		solution = np.append(solution, schr1c)
		solution = np.append(solution, schr2c)
		solution = np.append(solution, nhx)
		solution = np.append(solution, nhy)
		solution = np.append(solution, err_x)
		solution = np.append(solution, err_y)
		solution = np.append(solution, fobj)

		solution = np.reshape(solution, (1, len(solution)))
		solution = pd.DataFrame(solution)
		solution.to_csv('search.csv', mode='a', index=False, header=False)

	return fobj


########## Optimization ##################


import sys
import shutil


def _term_columns():
	return shutil.get_terminal_size((100, 20)).columns


def fitness_func(ga_instance, solution, solution_idx):
	fobj = objfunc(np.array(solution, dtype=float), print_table=False)

	gen_now = ga_instance.generations_completed + 1
	gen_total = ga_instance.num_generations

	try:
		pop_size = ga_instance.population.shape[0]
	except:
		pop_size = ga_instance.sol_per_pop

	try:
		current = int(solution_idx) + 1
	except:
		current = 0

	if pop_size > 0 and current > 0:
		percent = 100.0 * current / pop_size
		left = pop_size - current

		msg = (
			f"Generation {gen_now}/{gen_total} | "
			f"calculation {current}/{pop_size} | "
			f"left in generation: {left} | "
			f"{percent:6.2f}% | "
			f"f_obj = {fobj:.16e}"
		)

		width = _term_columns()
		print("\r" + msg[:width-1], end="", flush=True)

	return 1.0 / (1.0 + fobj)


def optimize(v0, num_generations=50, sol_per_pop=200, num_parents_mating=3,
             save_csv=True, csv_name="search_new.csv"):

	v1 = [np.array(v0, dtype=float)]
	nvars = len(v0)

	for i in range(500):
		v2 = np.zeros(nvars, dtype=float)

		if i == 0:
			for j in range(nvars):
				v2[j] = (1.0 + 0.0*np.random.normal()) * v0[j]
		elif i < 300:
			for j in range(nvars):
				v2[j] = (1.0 + 0.15*np.random.normal()) * v0[j]
		else:
			for j in range(nvars):
				v2[j] = (1.0 + 0.08*np.random.normal()) * v0[j]

		v1.append(v2)

	term_columns = _term_columns()

	print('=' * term_columns)
	print('INITIAL VALUES'.center(term_columns))
	print('=' * term_columns)

	fobj0 = objfunc(np.array(v0, dtype=float), print_table=True)

	print('=' * term_columns)
	print('Initial solution summary'.center(term_columns))
	print('=' * term_columns)
	print(f"{'Initial f_obj':<30} {fobj0:>25.16e}")
	print(f"{'Number of variables':<30} {nvars:>25}")
	print('=' * term_columns)
	print('Initial vector'.center(term_columns))
	print('=' * term_columns)
	print(np.array(v0, dtype=float))
	print('=' * term_columns)


	def on_generation(ga_instance):
		print()

		solution, solution_fitness, solution_idx = ga_instance.best_solution()
		fobj = objfunc(np.array(solution, dtype=float), print_table=True)

		term_columns = _term_columns()

		print('=' * term_columns)
		print(('GENERATION ' + str(ga_instance.generations_completed)).center(term_columns))
		print('=' * term_columns)

		print(f"{'Generation number':<30} {ga_instance.generations_completed:>25}")
		print(f"{'Best f_obj':<30} {fobj:>25.16e}")
		print(f"{'Best fitness':<30} {solution_fitness:>25.16e}")
		print(f"{'Best solution index':<30} {solution_idx:>25}")

		print('=' * term_columns)
		print('Genetic algorithm summary'.center(term_columns))
		print('=' * term_columns)

		try:
			pop_size = ga_instance.population.shape[0]
		except:
			pop_size = ga_instance.sol_per_pop

		print(f"{'Population size':<30} {pop_size:>25}")
		print(f"{'Number of genes':<30} {ga_instance.num_genes:>25}")
		print(f"{'Parents mating':<30} {ga_instance.num_parents_mating:>25}")
		print(f"{'Selection type':<30} {ga_instance.parent_selection_type:>25}")
		print(f"{'Crossover type':<30} {ga_instance.crossover_type:>25}")
		print(f"{'Mutation type':<30} {ga_instance.mutation_type:>25}")
		print(f"{'Mutation percent genes':<30} {ga_instance.mutation_percent_genes:>25}")

		print('=' * term_columns)
		print('Best solution vector'.center(term_columns))
		print('=' * term_columns)
		print(np.array(solution, dtype=float))
		print('=' * term_columns)

	ga_instance = pygad.GA(
		num_generations=num_generations,
		num_parents_mating=num_parents_mating,
		fitness_func=fitness_func,
		sol_per_pop=sol_per_pop,
		num_genes=nvars,
		initial_population=v1,
		init_range_low=-850,
		init_range_high=850,
		parent_selection_type="sss",
		crossover_type="single_point",
		mutation_type="random",
		mutation_percent_genes=20,
		on_generation=on_generation
	)

	ga_instance.run()

	print()

	solution, solution_fitness, solution_idx = ga_instance.best_solution()
	fobj = objfunc(np.array(solution, dtype=float), print_table=False)

	if save_csv:
		row = np.append(np.array(solution, dtype=float), fobj)
		pd.DataFrame(row.reshape(1, -1)).to_csv(csv_name, mode="a", index=False, header=False)

	term_columns = _term_columns()

	print('=' * term_columns)
	print('FINAL BEST SOLUTION'.center(term_columns))
	print('=' * term_columns)
	print(f"{'Best fitness':<30} {solution_fitness:>25.16e}")
	print(f"{'Best objective f_obj':<30} {fobj:>25.16e}")
	print(f"{'Best solution index':<30} {solution_idx:>25}")
	print('=' * term_columns)
	print(solution)
	print('=' * term_columns)

	return solution, fobj, ga_instance

########################################################
#EXECUTABLE FUNCTIONS
########################################################3

def main():
	global cell, linpart, ntel, ncel, term_size, energy, chrx, chry, schr1, schr2, whichprocess, varelem, varelemplus, cellfile, madxfile
	global stP, csftP, disptP, bxP
	global idx_to_vec, vec_to_idx, bracket_pairs, G, C, HatC, H_vec, M_basis, order, H_dict, Gnn
    #user setting
	m = 8
	d = 0
	n = 2  #Fixed for now!
	a_box = np.array([3e-2, 5e-3, 1.5e-3, 0.8e-3, 1e-3], dtype=float)
	H_pol = (
    sp.Rational(1, 2)*(px**2 + py**2)*(1 - delta + delta**2)
    - b1*x*delta
    + sp.Rational(1, 2)*b1**2*x**2
    + sp.Rational(1, 2)*b2*(x**2 - y**2)
    + sp.Rational(1, 3)*b3*(x**3 - 3*x*y**2)
    + sp.Rational(1, 4)*b4*(x**4 - 6*x**2*y**2 + y**4)
	)
	whichprocess='linear'
	opafile='opa_X_04-Feb-2025-10h53m_renglon_1.txt'
	#opafile='esrf.txt'
	cellfile='esrf.madx'
	madxfile='trackingone'
	chrx=0
	chry=0
	schr1='d2'
	schr2='d3'
	varelem=['qf1','qd2','qf3','qf1i','qd2i','qf3i','qf1d','qd2d','qf3d','d1','d4','d3d','d4d','o2','o3','o4','o3d','o4d']
	varelemplus=[]
	varelemplus=varelem.copy()
	varelemplus.append(schr1)
	varelemplus.append(schr2)
	varelemplus = ['k' + item for item in varelemplus]
	varelemplus.append('momentum')
	varelemplus.append('xamplitude')
	analysis='cell'
	linpart=0.01
	v0=[3.742153457895532e+00,-7.933030565890353e+00,9.383357074512245e+00,-1.961256869785227e+00,-7.622311324145445e+00,9.502339759326459e+00,1.528008695772636e+00,-6.400798167795065e+00,1.049016046874180e+01,3.148346595634115e+03*3,-4.621730496569649e+02*3,1.366797632793281e+02*3,-1.713849194855147e+02*3,-2.526041097299107e+12*1e-8,4.591682297815766e+07*1e-8,-3.136093044673090e+11*1e-8,-7.400633479124564e+10*1e-8,6.886521769314347e+11*1e-8]
	
	
    #loading
	ring,elements,energy=oparing(opafile)
	idx_to_vec, vec_to_idx, bracket_pairs, G, C, HatC, H_vec, M_basis, order, H_dict = load(m, d, H_pol, a_box, n=2)
	for i in ring:
		if i[0][0] == 'ring':
			ntel = len(i[1])
		if i[0][0] == analysis:
			ncel = len(i[1])

    #linear
	cell, csfunc, dispfunc, stP, csftP, disptP, bx, ax, gx, by, ay, gy, disp, dispd = build_linear_cell_and_functions(ring, elements, analysis, energy, linpart, whichprocess, plotlinear=0)
	bxP = bx
	#Si son globales aparentemente no tienes que recupeerarlas

	Gnn=G[15:,15:]

	global idx_x2, idx_xpx, idx_px2, idx_y2, idx_ypy, idx_py2, quad_size, nonquad_size
	idx_x2  = vec_to_idx[(0, 2, 0, 0, 0)]
	idx_xpx = vec_to_idx[(0, 1, 0, 1, 0)]
	idx_px2 = vec_to_idx[(0, 0, 0, 2, 0)]

	idx_y2  = vec_to_idx[(0, 0, 2, 0, 0)]
	idx_ypy = vec_to_idx[(0, 0, 1, 0, 1)]
	idx_py2 = vec_to_idx[(0, 0, 0, 0, 2)]

	quad_size = 15
	nonquad_size = len(idx_to_vec)-15

	fobj = objfunc(v0, print_table=True)
	quad_non_quad(bx, ax, gx, by, ay, gy)
	print("Fobj:")
	print(fobj)

	# quad_non_quad(bx, ax, gx, by, ay, gy)

	solution, fobj, ga_instance = optimize(v0, 5, 20)

	print("Final solution: ")
	print(solution)
	print("Final fobj:")
	print(fobj)
	#solution, fobj, ga_instance = optimize(v0,5,20)

	#print("Final solution: ")
	#print(solution)
	

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

    nhx, nhy, err_x, err_y= invariant(tnn, tnq, Sx, Sy, tol=1e-12)

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
	