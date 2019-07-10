#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright Serge Dmitrieff
# www.biophysics.fr

from numpy import *
import sys
import yaml
import pandas as pd
import time

__VERSION__ = "0.0.2"

"""
# SYNOPSIS
    Discrete reaction of actin monomers diffusing in a box and meeting the + tip

# DESCRIPTION
    The space is discretized as a 3D lattice and monomers move to a random adjacent time at each time step.
    Forbidden moves (collisions with a filament or the box boundary) are cancelled.

    Actin filaments can be modeled as :
    - Single filaments : orientation = [0,0]
    - Double filaments : orientation = [1,0] or [0,1]
    - Double filaments with a diffusion lattice half of the monomer size : orientation = [2,0] or [0,2]

# SYNTAX
    $ ./actin_polymerization.py [CONFIG_FILE] [OPTIONS]

    CONFIG_FILE is a suitable yaml-type config file, e.g.

        simulation:
            Tmax: 1000
            N_monomers: 100000
            box: [150,120,120]
            success_frac: 1.0

        filaments:
            f0:
                position: [60,19]
                orientation: [1,0]
                length: 50
            f1:
                position: [60,20]
                orientation: [1,0]
                length: 100

# OPTIONS
    See config file ; no other options available for now.

# EXAMPLES :
    $ ./actin_polymerization.py
    $ ./actin_polymerization.py config.cym

# OUTPUT :
    The output is a CSV file res.csv with the organization (with the aforementioned config file) :
    ,Tmax,N_monomers,box,success_frac,f0,f0.position,f0.orientation,f0.length,f1,f1.position,f1.orientation,f1.length
    ,1000,100000,"[150, 120, 120]",1.0,2847,"[60, 19]","[1, 0]",50,3688,"[60, 20]","[1, 0]",100

    In which the columns f0 and f1 yields the number of successful polymerization events for filaments f0 and f1 respectively

# @TODO :
    Implement change of length with polymerization (better simulation of diffusion limited growth)
    Implement an export for length(t)

"""


def version():
    return __VERSION__

class Integrator:
    """
        Class containing the Integrator

        Reads config file and creates a System.
    """

    # Constructor
    def __init__(self,*args,fname='config.yaml'):
        self.config=read_config(fname)
        self.system=System(*args,config=self.config)

    # Actual initializeation
    def initialize(self,*args,**kwargs):
        """ Initializes the integrator and the system """
        # First initializing the system
        self.system.initialize(*args,**kwargs)
        random.seed(int(time.time()))


        # Now we prepare the integrator using info from the config file
        config_sim=self.config['simulation']
        self.Tmax=config_sim['Tmax']
        self.N_filaments=self.system.N_filaments
        self.counts=zeros((self.N_filaments,),dtype= uint32)

    # Here we actually integrate the system with time
    def do_simulation(self,*args,**kwargs):
        for t in range(self.Tmax):
            self.counts+=self.system.make_step()

    # Saving result in csv format
    def save_results(self,*args,out="res.csv",**kwargs):
        # Results are stored in a dictionary
        results={}
        config_sim=self.config['simulation']

        # We record all the config values for the simulation
        for key in config_sim.keys():
            results[key]=config_sim[key]

        # We record all the config values for each filament
        for i,filament in enumerate(self.system.filaments):
            name=filament.name
            results[name]=self.counts[i]
            conf=filament.config
            for key in conf.keys():
                k='%s.%s' %(name,key)
                results[k]=conf[key]

        # We save it all in csv file format
        frame=pd.DataFrame(columns=results.keys())
        frame.loc['']=results
        frame.to_csv(out)
        return




class System:
    """
        Class containing the System

        Create all the array containing simulation state,
        and performs the time evolution step by step

        A system contains one or several Filaments
    """

    # Constructor
    def __init__(self,*args,config=None,**kwargs):
        if not config:
            raise ValueError('System cannot be initialized without a valid configuration')

        # Saving config in internal variables
        self.sim_config=config['simulation']
        self.fil_config=config['filaments']

        # Initializing boxes for convenience
        self.box=zeros((3,), dtype = int16)
        self.box_size=zeros((3,), dtype = int16)

        return

    # Loads the configuration to internal variables
    def load_config(self):
        conf=self.sim_config
        fils=self.fil_config
        # Loading from parameters
        self.N_monomers=int(conf['N_monomers'])
        self.filaments=[Filament(fils[key],name=key) for key in fils.keys()]
        self.N_filaments=len(self.filaments)
        self.box[:]=conf['box']
        self.box_size[:]=self.box[:]+1
        self.success_frac=conf['success_frac']
        self.positions=zeros((0,3),dtype= int16)


    def initialize(self,*args,**kwargs):
        # Actual system initialization
        self.load_config()

        # For further use
        Nm=self.N_monomers
        box=self.box

        # Initializing filaments
        for filament in self.filaments:
            filament.initialize(N_monomers=Nm,box=box)

        # box, but vectorized to the number of monomers
        self.big_box=full((Nm,1 ), 1, dtype=int16)*box

        # bool containers for convenience
        self.on=full((Nm, ), True, dtype=bool)
        self.out=full((Nm, 3), True, dtype=bool)
        self.bool_container=full((Nm, ), True, dtype=bool)

        #int container
        self.int_container=ones((Nm,),dtype=int16);

        # Computing random positions, setting displacement to 0
        self.confirmed=full((Nm, ), True, dtype=bool)
        self.displacement=zeros((Nm,3),dtype=int16)
        self.positions=self.random_position_inside(N=Nm)

        # setting counts to 0
        self.counts=zeros((self.N_filaments,),dtype= uint32)


    def random_position_inside(self,*args,N=1,**kwargs):
        # Returns an array of random positions inside the box !
        return hstack( (self.randcol_in(N=N,axis=0),self.randcol_in(N=N,axis=1),self.randcol_in(N=N,axis=2) ) )

    def randcol_in(self,N=1,axis=0):
        # Return a random of ints between 0 and box_size[axis]
        return random.randint(self.box_size[axis],size=(N,1),dtype=int16)

    def make_step(self):
        # step step step step
        self.compute_random_step()
        self.positions+=self.displacement

        # Checking who's out... And bringing them back in
        self.out[:,:]=self.positions<0
        self.positions[self.out]=0
        self.out[:,:]=self.positions>self.big_box
        self.positions[self.out]-=1

        # Now collision and reaction detection for each filament
        for i,filament in enumerate(self.filaments):
            # Reaction detection
            self.on[:]=((sum(self.positions==filament.big_pos[0],axis=1)==3)+(sum(self.positions==filament.big_pos[1],axis=1)==3))*self.confirmed
            N=sum(self.on)
            self.positions[self.on,:]=self.random_position_inside(N=N)
            self.counts[i]=N
            # Collision detection
            self.on[:]= sum([sum(self.positions[:,1:3]==fil,axis=1)==2 for fil in filament.big_fil ],axis=0)*(self.positions[:,0]<filament.length)
            self.positions[self.on,:]-=self.displacement[self.on]

        return self.counts


    def compute_random_step(self):
        # Here we compute the random variables
        self.displacement.fill(0)
        self.int_container[:]=random.randint(3,size=(self.N_monomers,),dtype=int16)
        for i in range(3):
            j=int16(i)
            self.bool_container[:]=(self.int_container[:]==j)
            self.displacement[self.bool_container,i]=2*random.randint(2,size=(sum(self.bool_container),),dtype=int16)-1
        if self.success_frac<1:
            self.confirmed[:]=random.binomial(1,self.success_frac,self.N_monomers)


class Filament:
    """
        Class containing the Filament

        Reads config file and creates a filament.
    """

    # Constructor
    def __init__(self,config,name='filament',**kwarg):
        self.config=config
        self.position=array(config['position'],dtype=int16)
        self.orientation=array(config['orientation'],dtype=int16)
        self.length=int16(config['length'])
        self.name=name

    # actual initialization
    def initialize(self,N_monomers=0,box=array([0,0,0]),**kwargs):
        """ Here we initialize the filament coordinates

            big_fil is a list of Nm x 2 arrays,
                made to rapidly compare monomer coordinates to filament coords in the YZ plane
                in order to detect collisions

            big_fil is a list of Nm x 3 arrays,
                made to rapidly compare monomer coordinates the + tip position
                in order to detect collisions
        """

        Nm=N_monomers
        # Coords contain the coordinates of 'protofilaments',
        #   e.g. contiguous lines occupied by filaments
        coords=self.make_coordinates()
        # position in the YZ plane
        position=self.position
        # Orientation in the YZ plane
        dir=self.orientation
        length=self.length
        # XYZ position of the tip
        posis=[hstack((array([length],dtype=int16),position)),hstack((array([length],dtype=int16),position+dir))]
        # big_pos is to dectect collisions along the filament
        self.big_pos=[full((Nm, 1), 1, dtype=int16)*pos for pos in posis  ]
        # big_pos is to dectect collisions at the tip and potential polymerization
        self.YZ=[array(position+coord,dtype=int16) for coord in coords]
        self.big_fil=[full((Nm, 1), 1, dtype=int16)*yz  for yz  in self.YZ]

    # There we compute the position of the lines (along x) occupied by filaments
    def make_coordinates(self):
        coords=[]
        dir=self.orientation
        if abs(sum(dir))==2:
            # DOUBLE RESOLUTION !
            if abs(dir[0])>0:
                for x in range(-1,4):
                    for y in range(-1,2):
                        coords.append([x,y])
            else:
                for x in range(-1,2):
                    for y in range(-1,4):
                        coords.append([x,y])
        elif abs(sum(dir))==1:
            # Simple resolution, 1 filament = 2 protofilaments
            coords.append([0,0])
            coords.append(dir)
        else:
            # Most likely single filaments
            coords.append(dir)
        return array(coords,dtype=int16)


def read_config(fname):
    # Reads a config file to a dictionary
    #   a wrapper to yaml.load
    file=open(fname,'r')
    config=yaml.load(file, Loader=yaml.SafeLoader)
    file.close()
    return config

if __name__ == "__main__":
    # The program as called from the command line
    #   Here we just gather arguments to create and run the Integrator

    # the config file, if any, should be the first argument
    nargs=len(sys.argv)
    args=[]
    if nargs<2:
        # if no config file name was given
        fname="config.cym"
    else:
        # a config file name was given !
        fname=sys.argv[1]
        if nargs>2:
            # if there are more arguments (there shouldn't for now)
            args=sys.argv[2:]

    # Now we create the integrator
    integ=Integrator(*args,fname=fname)
    # We initialize and run the simulation
    integ.initialize()
    integ.do_simulation()
    # Also, saving is good
    integ.save_results()
