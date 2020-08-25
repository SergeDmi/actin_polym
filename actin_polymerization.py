#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright Serge Dmitrieff
# www.biophysics.fr

#from numpy import *
import numpy as np
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

    The position of a filament is that of its TOP LEFT corner.

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

    The position of the filament is the position of the top left corner in the yz plane

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
    Implement an export for length(t) (or at least length(Tend) )

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
        self.system=System(*args,integrator=self,config=self.config)

    # Actual initializeation
    def initialize(self,*args,**kwargs):
        """ Initializes the integrator and the system """
        # First initializing the system
        self.system.initialize(*args,**kwargs)
        np.random.seed(int(time.time()))
        self.current_t=0

        # Now we prepare the integrator using info from the config file
        config_sim=self.config['simulation']
        self.Tmax=config_sim['Tmax']
        self.N_filaments=self.system.N_filaments
        self.counts=np.zeros((self.N_filaments,),dtype= np.uint32)
        self.make_data_frame()

    # Here we actually integrate the system with time
    def do_simulation(self,*args,**kwargs):
        for t in range(self.Tmax):
            self.current_t=t
            self.counts+=self.system.make_step()

    def make_results(self):
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
                if key=="length":
                    results[k]=filament.get_length_from_config()
                    results['%s.final_%s' %(name,key)]="%s" %filament.length
                else:
                    results[k]=conf[key]
        return results

    # Creates a data frame to save results in
    def make_data_frame(self):
        results=self.make_results()
        self.frame=pd.DataFrame(columns=results.keys())

    # Saving result in csv format
    def save_current_state(self,*args,out="res.csv",**kwargs):
        # Results are stored in a dictionary
        results=self.make_results()
        # We save it all in csv file format
        self.frame.loc[self.current_t]=results

        return


    def export_results(self,*args,out="res.csv",**kwargs):
        # Results are stored in a dictionary
        self.frame.to_csv(out)
        return

class System:
    """
        Class containing the System

        Create all the array containing simulation state,
        and performs the time evolution step by step

        A system contains one or several Filaments
    """

    # Constructor
    def __init__(self,*args,config=None,integrator=None,**kwargs):
        if not config:
            raise ValueError('System cannot be initialized without a valid configuration')

        self.integrator=integrator
        # Saving config in internal variables
        self.sim_config=config['simulation']
        self.fil_config=config['filaments']

        # Initializing boxes for convenience
        self.box=np.zeros((3,), dtype = np.int16)
        self.box_size=np.zeros((3,), dtype = np.int16)

        self.dynamic_filaments=0
        return

    # Loads the configuration to internal variables
    def load_config(self):
        conf=self.sim_config
        fils=self.fil_config
        # Loading from parameters
        self.N_monomers=int(conf['N_monomers'])
        self.filaments=[Filament(fils[key],system=self,name=key) for key in fils.keys()]
        self.N_filaments=len(self.filaments)
        self.box[:]=conf['box']
        self.box_size[:]=self.box[:]+1
        self.success_frac=conf['success_frac']
        self.positions=np.zeros((0,3),dtype= np.int16)
        self.arange=np.arange(self.N_monomers)
        try:
            self.dynamic_filaments=conf['growing']
        except:
            self.dynamic_filaments=0


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
        self.big_box=np.full((Nm,1 ), 1, dtype=np.int16)*box

        # bool containers for convenience
        self.on=np.full((Nm, ), True, dtype=np.bool)
        self.out=np.full((Nm, 3), True, dtype=np.bool)
        self.bool_container=np.full((Nm, ), True, dtype=np.bool)

        #int container
        self.int_container=np.ones((Nm,),dtype=np.int16);

        # Computing random positions, setting displacement to 0
        self.confirmed=np.full((Nm, ), True, dtype=np.bool)
        self.displacement=np.zeros((Nm,3),dtype=np.int16)
        self.positions=self.random_position_inside(N=Nm)

        # setting counts to 0
        self.counts=np.zeros((self.N_filaments,),dtype= np.uint32)

    def random_position_inside(self,*args,N=1,**kwargs):
        # Returns an array of random positions inside the box !
        #   Now avoids the filaments
        # Warning : most of the times, this is called with very very few N :
        # probably not performance limiting
        done=0
        left=self.arange[0:N]
        pos=np.zeros((N,3),dtype=np.int16)
        to_do=N

        while to_do:
            # We assume no monomer is on a filament
            on=np.full((to_do, ), False, dtype=np.bool)
            # We make random positions for monomers
            pos[left,:]=np.hstack( (self.randcol_in(N=to_do,axis=0),self.randcol_in(N=to_do,axis=1),self.randcol_in(N=to_do,axis=2) ) )
            # we check if any monomer is on a filament
            for i,filament in enumerate(self.filaments):
                on[:]+=np.array( np.sum([ ( np.sum(pos[left,1:3]==sub_fil[left,:],axis=1)==2 )*( pos[left,0] < filament.sub_lengths[j] ) for j,sub_fil in enumerate(filament.big_fil) ],axis=0) , dtype=np.bool)
            # we select monomers that are actually on a filament
            left=left[on]
            to_do=len(left)

        return pos


    def randcol_in(self,N=1,axis=0):
        # Return a random of ints between 0 and box_size[axis]
        return np.random.randint(self.box_size[axis],size=(N,1),dtype=np.int16)

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

            # Collision detection
            self.on[:]= np.sum([ ( np.sum(self.positions[:,1:3]==sub_fil,axis=1)==2 )*( self.positions[:,0] < filament.sub_lengths[j] ) for j,sub_fil in enumerate(filament.big_fil) ],axis=0)
            #self.on[:]= sum([(sum(self.positions[:,1:3]==fil,axis=1)==2)*(self.positions[:,0]<filament.length_fil[j]) for j,fil in enumerate(filament.big_fil)],axis=0)
            # if collision, we cancel diffusion
            self.positions[self.on,:]-=self.displacement[self.on]



            # Reaction detection
            active=filament.reactive
            self.on[:]=np.sum(self.positions==filament.big_pos[active],axis=1)==3
            #b=filament.big_pos[active]
            #print(b[0,:])
            N=np.sum(self.on)

            # We check if actual reaction, or fail
            if self.success_frac<1:
                self.on[:]*=np.random.binomial(1,self.success_frac,N)
                N=np.sum(self.on)

            self.counts[i]=N

            if N>0:
                if self.dynamic_filaments:
                    # If filaments are dynamic, only a single polymerization event can happen at a given time
                    N=1
                    # If filaments are dynamic we make them grow !
                    filament.add_length(N)
                    ix=min(self.arange[self.on])
                    # monomers that reacted are injected back in the solution
                    self.positions[ix,:]=self.random_position_inside(N=N)
                    # We check if reaction didn't cause new collisions. if so, we push monomers around
                    self.on[:]= np.sum([ ( np.sum(self.positions[:,1:3]==sub_fil,axis=1)==2 )*( self.positions[:,0] < filament.sub_lengths[j] ) for j,sub_fil in enumerate(filament.big_fil) ],axis=0)
                    self.positions[self.on,1:3]+=filament.brick



        return self.counts

    def compute_random_step(self):
        # Here we compute the random variables
        self.displacement.fill(0)
        self.int_container[:]=np.random.randint(3,size=(self.N_monomers,),dtype=np.int16)
        for i in range(3):
            j=np.int16(i)
            self.bool_container[:]=(self.int_container[:]==j)
            self.displacement[self.bool_container,i]=2*np.random.randint(2,size=(np.sum(self.bool_container),),dtype=np.int16)-1
        #if self.success_frac<1:
        #    self.confirmed[:]=random.binomial(1,self.success_frac,self.N_monomers)

class Filament:
    """
        Class containing the Filaments

        a filament contains a set of "sub-filaments" != proto-filaments
        the set of sub-filaments correspond to space around the filament where monomers cannot be

        1 Proto-filament :

        X

        2 Proto-filaments, grid size=monomer size : sub-filament = proto-filament

        X     or      X X
        X

        2 Proto-filaments, grid size=monomer size/2
            -> 15 sub-filaments

        sub-filaments can be visualized as :
        _ _ _         _ _ _ _ _
        _ 0 X    or   _ 0 X 0 X
        _ X X         _ X X X X
        _ 0 X
        _ X X

        Where "X", '_' and "0" denotes the locations of sub-filaments (...)
        "0" denotes the subfilament to which monomer can react
        "X" and "0" is where actin actually is
        "_","X","0" are forbidden positions for monomers
        A monomer looks something like this :
        _ _ _
        _ 0 X
        _ X X


    """

    # Constructor
    def __init__(self,config,name='filament',system=None,**kwarg):
        self.system=system
        self.config=config
        self.position=np.array(config['position'],dtype=np.int16)
        self.orientation=np.array(config['orientation'],dtype=np.int16)
        self.brick=max(abs(self.orientation))

        self.length=self.get_length_from_config()
        self.name=name
        self.reactive=self.get_reactive()

    # actual initialization
    def initialize(self,N_monomers=0,box=np.array([0,0,0]),**kwargs):
        """
            Here we initialize the filament coordinates

            big_fil is a list of Nm x 2 arrays,
                made to rapidly compare monomer coordinates to filament coords in the YZ plane
                in order to detect collisions

            big_pos is a list of Nm x 3 arrays,
                made to rapidly compare monomer coordinates the + tip position
                in order to detect collisions
        """

        Nm=N_monomers
        # coords contain the coordinates of 'protofilaments',
        #   e.g. contiguous lines occupied by filaments
        # dir contains the orientation in the YZ plane
        coords,dir=self.make_coordinates()
        # position in the YZ plane
        #print("coords:")
        #print(coords)
        position=self.position
        # length of filament
        length=self.length

        # XYZ position of the tip
        posis=[np.hstack((np.array([length[0]],dtype=np.int16),position)),np.hstack((np.array([length[1]],dtype=np.int16),position+dir))]

        # big_pos is to dectect collisions along the filament
        self.big_pos=[np.full((Nm, 1), 1, dtype=np.int16)*pos for pos in posis  ]
        # YZ containts the XY positions of single mini filaments
        self.YZ=[np.array(position+coord,dtype=np.int16) for coord in coords]


        self.big_fil=[np.full((Nm, 1), 1, dtype=np.int16)*yz  for yz  in self.YZ]

        # initiation with 0
        self.sub_lengths=0.0*coords[:,0]
        # fill in the sub lengths
        self.update_lengths()

        if self.system.dynamic_filaments:
            self.system.integrator.save_current_state()

    # There we compute the position of the lines (along x) occupied by filaments
    def make_coordinates(self):
        """
            Here we make the coordinates of all sub-filaments of a given filament
        """
        coords=[]
        dir=self.orientation
        if abs(np.sum(dir))==2:
            # DOUBLE RESOLUTION !
            # Now assuming filament position is top left
            if abs(dir[0])>0:
                for x in range(-1,4):
                    for y in range(-1,2):
                        coords.append([x,y])
                dir=np.array([2,0])
            else:
                for y in range(-3,2):
                    for x in range(-1,2):
                        coords.append([x,y])
                dir=np.array([0,-2])
        elif abs(np.sum(dir))==1:
            # Simple resolution, 1 filament = 2 protofilaments
            if abs(dir[1])==1:
                dir=np.array([0,-1])
            else:
                dir=np.array([1,0])
            coords.append([0,0])
            coords.append(dir)
        else:
            # Most likely single filaments
            coords.append(dir)
        return np.array(coords,dtype=np.int16),dir

    # finds out which end is reactive (the shortest is)
    def get_reactive(self):
        return np.int16(self.length[1]<self.length[0])

    # Understanding the length as an array
    def get_length_from_config(self):
        config=self.config
        if self.brick:
            try:
                length=np.array([np.int16(config['length']),np.int16(config['length'])],dtype=np.int16)
            except:
                try:
                    length=np.array(config['length'],dtype=np.int16)
                except:
                    raise ValueError('Could not understand length from %s' % config['length'])
            if length[0]==length[1]:
                length[0]+=1
        else:
            try:
                length=np.int16(config['length'])
            except:
                raise ValueError('Could not understand length from %s' % config['length'])
        return length

    def add_length(self,N):
        self.length[self.reactive]+=N*self.brick
        self.update_lengths()
        self.system.integrator.save_current_state()
        #print(self.length)

    # update lengths of subfilaments
    def update_lengths(self):
        """
            Here we update the length of sub filaments

                Numerotation of sub filaments :

                (vertical)          (horizontal)

                12 13 14            2 5 8 11 14
                9  10 11            1 4 7 10 13
                6  7  8             0 3 6 9  12
                3  4  5
                0  1  2

                So reaction can happen at sub-filaments 4 or 10.
        """

        self.reactive=self.get_reactive()

        # The sub-filaments have different lengths
        if self.brick<2:
            self.sub_lengths=self.length
        else:
            # here it's important to remember that one protofilament is
            if self.reactive==0:
                self.sub_lengths[0:6]=self.length[0]
                self.sub_lengths[6:]=self.length[1]
            else:
                self.sub_lengths[0:9]=self.length[0]
                self.sub_lengths[9:] =self.length[1]

        for i,pos in enumerate(self.big_pos):
            pos[:,0]=self.length[i]


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
    integ.save_current_state()
    integ.export_results()
