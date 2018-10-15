# coding: utf-8

import os
import shutil
from shutil import rmtree
from monty.shutil import decompress_dir, decompress_file

from pymatgen.io.vasp import Vasprun

from pawpyseed.core.wavefunction import Wavefunction
from pawpyseed.core.projector import Projector

class PathHolder():
	def __init__(self, path):
		self.launch_dir = path

class DummyFirework():
    def __init__(self, path):
        self.launches = [PathHolder(path)]
        self.name = path
        self.fw_id = path

class DefectWorkflowWavefunctionHandle(object):
    """
    This class is made to run Kyle's
    wavefunction parsing code...

    bulk_fw_sets is a dict that has bulk fireworks with supercell size as key
    wf_job_run is a string which says if it is "normal" (=gga) or "scan" or "hse"
    """
    def __init__(self, bulk_fw_sets, dwo = None):
        self.bulk_fw_sets = bulk_fw_sets
        self.dwo = dwo

    def _setup_file_for_parsing(self, path):
        #make a "kyle_file" for parsing out relevant data
        #returns True if file setup correctly
        files_to_copy = ['CONTCAR','OUTCAR','POTCAR',
                         'WAVECAR','vasprun.xml']

        kyle_path = os.path.join( path, 'kyle_file')
        if os.path.exists(kyle_path):
            print('KYLE FILE ALREADY! Removing and rebuilding...')
            rmtree(kyle_path)

        #make kyle file
        os.makedirs(kyle_path)

        #copy relevant files to kyle path
        for file in files_to_copy:
            filepat = os.path.join( path, file +'.relax2.gz')
            if not os.path.exists( filepat):
                filepat = os.path.join( path, file +'.relax1.gz')
            if not os.path.exists( filepat):
                filepat = os.path.join( path, file +'.gz')
            if not os.path.exists( filepat):
                filepat = os.path.join( path, file)
            if not os.path.exists( filepat):
                print('Could not find {}! Skipping this defect...'.format(file))
                return False

            #COPY file to kyle_file
            if '.gz' in filepat:
                shutil.copy(filepat, os.path.join( kyle_path, file+'.gz'))
            else:
                shutil.copy(filepat, os.path.join( kyle_path, file))

        decompress_dir( kyle_path)
        return True

    def _get_vbm_band_dict(self, path, vbm):
        # ASSUMING KYLE PATH ALREADY SET UP, RETURNS A BAND_DICT to be used by run_pawpy
        # (band numbers are keys, contains maxmin window of band energies (along with occupation),
        #  and stores percentage of band character after pawpy is run...)
        vr = Vasprun(os.path.join(path, 'kyle_file', 'vasprun.xml'))
        max_num = 0
        band_dict = {bandindex: {'max_eigen': [-10000., 0.], 'min_eigen': [10000., 0.],
                                'VB_projection': None, 'CB_projection': None}
                                for bandindex in range(len(list(vr.eigenvalues.values())[0][0]))}

        for spin, spinset in vr.eigenvalues.items():
            for kptset in spinset:
                for bandnum, eigenset in enumerate(kptset):
                    if eigenset[1] and (eigenset[0] <= vbm) and (bandnum > max_num):
                        max_num = bandnum
                    #see if this is lowest eigenvalue for this band so far
                    if eigenset[0] < band_dict[bandnum]['min_eigen'][0]:
                        band_dict[bandnum]['min_eigen'] = eigenset[:]
                    #see if this is highest eigenvalue for this band so far
                    if eigenset[0] > band_dict[bandnum]['max_eigen'][0]:
                        band_dict[bandnum]['max_eigen'] = eigenset[:]

        trim_band_dict = {band_index: band_dict[band_index].copy() for band_index in range(max_num-20, max_num+21)}

        return trim_band_dict, vr.is_spin


    def run_pawpy(self):
        """
        Container for running pawpyseed on all defects in this workflow
        """

        vbm = None
        bulk_dirs, wf_dirs = [], []
        bulk_sizes, wf_sizes = [], []
        for sc_size, bulk_fw in self.bulk_fw_sets.items():
            launch_dir = bulk_fw.launches[-1].launch_dir
            bulk_sizes.append(sc_size)
            bulk_dirs.append(launch_dir)
            if not vbm:
                # need to check different filenames
                vr = Vasprun(os.path.join(launch_dir, 'vasprun.xml'))
                vbm = vr.eigenvalue_band_properties[2]
                print('\twill use vbm value of ',vbm)
        for sc_size, size_set in self.dwo.defect_fw_sets.items():
            for fw in size_set:
                wf_sizes.append(sc_size)
                wf_dirs.append(fw.launches[-1].launch_dir)
        projector_list, bases = Projector.setup_bases(bulk_dirs, wf_dirs, True)
        num_proj_els = bases[0].num_proj_els
        store_all_data = {}
        basis_sets = {}
        for i, sc_size in enumerate(bulk_sizes):
            basis_sets[sc_size] = bases[i]

        #for each defect, in a given supercell size
        for sc_size, size_set in self.dwo.defect_fw_sets.items():
            for fw in size_set:
                try:
                    launch_dir = fw.launches[-1].launch_dir
                    print('start parsing of {} wavefunction'.format(fw.name))

                    # setup file path,
                    print('\tsetting up files')
                    self._setup_file_for_parsing( launch_dir)

                    # find band number which is at VBM and iniitalize band_dict and find spin polarization
                    print('\tinitializing band_dict')
                    band_dict, spinpol = self._get_vbm_band_dict( launch_dir, vbm)

                    #setup defect wavefunction
                    print('\tmerging wf from dir')
                    wf = Wavefunction.from_atomate_directory( launch_dir, setup_projectors=False)

                    # loop over band projections around band edge and store results
                    print('\tperforming projections')
                    pr = Projector(wf, basis_sets[sc_size],
                        projector_list = projector_list, unsym_wf = True)
                    for bandnum in band_dict.keys():
                        v,c = pr.proportion_conduction( bandnum, spinpol=spinpol)
                        band_dict[bandnum]['VB_projection'] = v[:]
                        band_dict[bandnum]['CB_projection'] = c[:]
                        print(bandnum, band_dict[bandnum])

                    #then tear down file path
                    print('\ttear down files')
                    rmtree(os.path.join( launch_dir, 'kyle_file'))

                    #release defect memory
                    pr.wf.free_all()

                    store_all_data[fw.fw_id] = band_dict
                except Exception as e:
                    print("___&*$#&(*#@&$)(*&@#)($----\n--> ERROR OCCURED. "
                          "Skipping this defect.\n-------------^#$^*&^#$&*^#@^$#-------------")
                    print(repr(e))

        #now release all bulk basis memory
        for bulk_basis in basis_sets.values():
            bulk_basis.free_all()

        Projector.free_projector_list(projector_list, num_proj_els)

        return store_all_data
