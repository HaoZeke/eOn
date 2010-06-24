import os
import shutil
import logging
logger = logging.getLogger('communicator')

from time import sleep, time
import subprocess
import commands
import tarfile

class NotImplementedError(Exception):
    pass

class CommunicatorError(Exception):
    pass

class Communicator:
    def __init__(self, scratchpath, bundle_size=1):
        if not os.path.isdir(scratchpath):
            #should probably log this event
            os.makedirs(scratchpath)
        self.scratchpath = scratchpath
        self.bundle_size = bundle_size

    def submit_searches(self, jobpaths):
        '''Throws CommunicatorError if fails.'''
        raise NotImplementedError()

    def get_results(self, resultspath):
        '''Returns a list of directories containing the results.'''
        raise NotImplementedError()

    def get_queue_size(self):
        '''Returns the number of items waiting to run in the queue.'''
        raise NotImplementedError()

    def cancel_state(self, statenumber):
        '''Returns the number of workunits that were canceled.'''
        raise NotImplementedError()

    def get_bundle_size(self, result_dat_path):
        try:
            f = open(result_dat_path)
        except:
            logger.warning("results.dat missing. Should be at %s.", result_dat_path)
            return 0
        bundle_size = 0
        for line in f:
            fields = line.split()
            if fields[1] == "termination_reason":
                bundle_size += 1
        if bundle_size == 0:
            logger.warning("termination_reason missing for result file %s", result_dat_path)
        return bundle_size

    def unbundle(self, resultpath):
        '''This method unbundles multiple searches into multiple single 
           searchs so the akmc script can process them.'''

        # These are the files in the result directory that we keep.
        file_list = [ "mode.dat", "product.con", "reactant.con",
                      "results.dat", "saddle.con"]

        jobpaths = [ os.path.join(resultpath,d) for d in os.listdir(resultpath) 
                    if os.path.isdir(os.path.join(resultpath,d)) ]

        result_dirs = []
        for jobpath in jobpaths:
            # Need to figure out how many searches were bundled together
            # and then create the new job directories with the split files.
            bundle_size = self.get_bundle_size(os.path.join(jobpath, "results.dat"))
            # Get the number at the end of the jobpath. Looks like path/to/job/#_#
            # and we want the second #.
            basename, dirname = os.path.split(jobpath)
            state, uid = dirname.split('_')
            state = int(state)
            uid = int(uid)
            jobpath_template = os.path.join(basename, "%i_%i")
            fdata = {}
            for i in range(bundle_size):
                jobpath_dest = jobpath_template % (state, uid+i)
                # Read in the files we care about into memory so we can split them
                # and then delete the files.
                for fname in file_list:
                    fpath = os.path.join(jobpath_dest, fname)
                    if i == 0:
                        f = open(fpath)
                        fdata[fname] = f.readlines()
                        f.close()
                        os.unlink(fpath)
                    else:
                        if not os.path.isdir(jobpath_dest):
                            os.makedirs(jobpath_dest)

                    #check to see that len(fdata[fname]%bundle_size == 0
                    offset = len(fdata[fname])/bundle_size
                    f = open(fpath, "w")
                    fsplitdata = ''.join(fdata[fname][i*offset:(i+1)*offset])
                    f.write(fsplitdata)
                    f.close()

                result_dirs.append(os.path.split(jobpath_dest)[1])
        return result_dirs

    def make_bundles(self, jobpaths):
        '''This method is a generator that bundles together multiple searches into a single job.
           Example usage:
               for jobpath in self.make_bundles(jobpaths):
                   do_stuff()'''
        # Split jobpaths in to lists of size self.bundle_size.
        chunks = [ jobpaths[i:i+self.bundle_size] for i in range(0, len(jobpaths), self.bundle_size) ]
        for chunk in chunks:
            # Open the first jobpath's displacement and mode files.
            dp_concat = open(os.path.join(chunk[0],"displacement_passed.con"), "a")
            mp_concat = open(os.path.join(chunk[0],"mode_passed.dat"), "a")

            # Concatenate all of the displacement and modes together.
            for i in range(1,len(chunk)):
                dp_join = open(os.path.join(jobpaths[i], "displacement_passed.con"))
                mp_join = open(os.path.join(jobpaths[i], "mode_passed.dat"))
                
                dp_concat.write(dp_join.read())
                mp_concat.write(mp_join.read())

                dp_join.close()
                mp_join.close()

            dp_concat.close()
            mp_concat.close()

            # Returns the jobpath to the new bigger workunit.
            yield chunk[0]

class BOINC(Communicator):
    def __init__(self, scratchpath, boinc_project_dir, wu_template, 
            result_template, appname, boinc_results_path, bundle_size):
        '''This constructor modifies sys.path to include the BOINC python modules. 
        It then tries to connect to the BOINC mysql database raising exceptions if there are problems 
        connecting. It also creates a file named uniqueid in the scratchpath to identify BOINC jobs as
        belonging to this akmc run if it doesn't exist. It then reads in the uniqueid file and 
        stores that as an integer in self.uniqueid.'''
        
        Communicator.__init__(self, scratchpath, bundle_size)
        self.wu_template = wu_template
        self.result_template = result_template
        self.appname = appname
        self.boinc_project_dir = boinc_project_dir
        self.boinc_results_path = boinc_results_path
        os.environ['BOINC_PROJECT_DIR'] = self.boinc_project_dir
        import sys
        sys.path.insert(0, os.path.join(self.boinc_project_dir, 'py'))
        try:
            import Boinc.database
            import Boinc.db_base
            import Boinc.boinc_db
        except ImportError:
            raise CommunicatorError("The Boinc python module could not be imported.\nPerhaps "
                "the boinc project path is set incorrectly?")

        self.database = Boinc.database 
        self.boinc_db_constants = Boinc.boinc_db
        self.db_base = Boinc.db_base

        try:
            self.database.connect_default_config()
        except:
            # XXX: This error handling is maybe a little ugly, but provides all the information
            # that you would want to know. The exception that connect_default_config() throws
            # is not helpful. It often will just say that it couldn't parse a xml file when it
            # really means it can't find the project's config.xml file.
            import traceback
            # print the traceback from connect_default_config
            traceback.print_exc()
            # raise a nice human readable error
            raise CommunicatorError("Couldn't connect to the BOINC database.")

        self.dbconnection = self.db_base.dbconnection
        self.cursor = self.dbconnection.cursor()

        #generate our unique id if it doesn't already exist.
        uniqueid_path = os.path.join(self.scratchpath, "uniqueid")

        if not os.path.isfile(uniqueid_path):
            f = open(uniqueid_path, 'w')
            import random 
            uid = random.randint(0, 2147483647)
            f.write("%s\n" % uid)
            f.close()
            logger.debug("wrote new unique id %i to %s" % (uid, uniqueid_path))
        try:
            f = open(uniqueid_path)
        except IOError:
            raise CommunicatorError("Unable to open the uniqueid file: %s" % uniqueid_path)

        try:
            self.uniqueid = int(f.read().strip())
            logger.debug("read in unique id %i from %s" % (self.uniqueid, 
                uniqueid_path))
        except ValueError:
            raise CommunicatorError("Trouble converting uniqueid value in %s to integer" % uniqueid_path)

    def get_queue_size(self):
        server_state = self.boinc_db_constants.RESULT_SERVER_STATE_UNSENT
        query = 'select count(*) from result where batch=%i and server_state=%i'
        query = query % (self.uniqueid, server_state)
        self.cursor.execute(query)
        row = self.cursor.fetchone()
        number_unsent = row['count(*)']

        return number_unsent

    def cancel_state(self, statenumber):
        # XXX: This function might be too expensive. Probably needs to be
        #      profiled later. It has to get all of the result rows that correspond
        #      to a unique id find which state they correspond to by parsing the
        #      wu name and then update the rows that correspond to statenumber.

        state_unsent = self.boinc_db_constants.RESULT_SERVER_STATE_UNSENT
        q1 = "select id,workunitid,name from result where batch=%i and server_state=%i" 
        q1 = q1 % (self.uniqueid, state_unsent)

        self.cursor.execute(q1)
        
        if self.cursor.rowcount == 0:
            return 0

        result_ids = []
        workunit_ids = []
        while 1:
            row = self.cursor.fetchone()
            if row == None:
                break
            
            result_id = row['id']
            workunit_id = row['workunitid']
            name = row['name']
            result_statenumber = name.split('_')[1]

            if statenumber == int(result_statenumber):
                result_ids.append(str(result_id))
                workunit_ids.append(str(workunit_id))

        resultid_string = '('+','.join(result_ids) + ')'
        state_over = self.boinc_db_constants.RESULT_SERVER_STATE_OVER
        outcome_not_needed = self.boinc_db_constants.RESULT_OUTCOME_DIDNT_NEED
        error_mask_cancelled = self.boinc_db_constants.WU_ERROR_CANCELLED
        q1 = "update result set server_state=%i, outcome=%i where id in %s" 
        q1 = q1 % (state_over, outcome_not_needed, resultid_string)
        self.cursor.execute(q1)

        workunitid_string = '('+','.join(workunit_ids) + ')'
        q2 = "update workunit set error_mask=%i, transition_time=%i where id in %s" 
        q2 = q2 % (error_mask_cancelled, int(time()), workunitid_string)
        self.cursor.execute(q2)
        num_cancelled_wu = self.cursor.rowcount

        self.db_base.dbconnection.commit()

        return num_cancelled_wu

    def dir_hier_path(self, filename):
        cmd = os.path.join(self.boinc_project_dir,"bin","dir_hier_path")
        path = commands.getoutput("%s %s" % (cmd, filename))
        return path

    def submit_searches(self, jobpaths):
        '''Runs the BOINC command create_work on all the jobs.'''
        from threading import Thread
        thread_list = []
        for jobpath in self.make_bundles(jobpaths):
            wu_name = "%i_%s" % (self.uniqueid, os.path.split(jobpath)[1])
            t = Thread(target=self.create_work, args=(jobpath, wu_name))
            t.start()
            thread_list.append(t)

            if len(thread_list) == 10:
                for t in thread_list:
                    t.join()
                thread_list = []

    def create_work(self, jobpath, wu_name):
        create_wu_cmd = os.path.join('bin', 'create_work')
        rp_path = self.dir_hier_path('reactant_passed_%s.con' % wu_name).strip()
        pp_path = self.dir_hier_path('parameters_passed_%s.dat' % wu_name).strip()
        dp_path = self.dir_hier_path('displacement_passed_%s.con' % wu_name).strip()
        mp_path = self.dir_hier_path('mode_passed_%s.dat' % wu_name).strip()

        shutil.move(os.path.join(jobpath, 'reactant_passed.con'), rp_path)
        shutil.move(os.path.join(jobpath, 'parameters_passed.dat'), pp_path)
        shutil.move(os.path.join(jobpath, 'displacement_passed.con'), dp_path)
        shutil.move(os.path.join(jobpath, 'mode_passed.dat'), mp_path)

        arglist = [create_wu_cmd]
        arglist.append("-appname")
        arglist.append(self.appname)
        arglist.append("-wu_name")
        arglist.append(wu_name)
        arglist.append("-wu_template")
        arglist.append(self.wu_template)
        arglist.append("-result_template")
        arglist.append(self.result_template)
        arglist.append("-batch")
        arglist.append(str(self.uniqueid))

        # XXX: This assumes an order of the input files. We should read the
        #      input template xml file to discover this order. Too bad the
        #      job templates aren't valid XML files. This means that a
        #      custom template file reader has to be written to read these
        #      files.
        arglist.append(os.path.split(rp_path)[1])
        arglist.append(os.path.split(pp_path)[1])
        arglist.append(os.path.split(dp_path)[1])
        arglist.append(os.path.split(mp_path)[1])

        #logger.debug("submited wu %s" % wu_name)
        p = subprocess.Popen(arglist, cwd=self.boinc_project_dir,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        #TODO: Benchmark this method. I previously found that spawning
        #      a bunch of create_work commands and waiting for all of them
        #      to finish in batches was faster. No premature optimization 
        #      here.
        retval = p.wait()

        if retval != 0:
            stdout, stderr = p.communicate()
            errstr = "Problem submitting the BOINC workunit:\nstdout: %s\nstderr %s" % \
                        (stdout, stderr)
            raise CommunicatorError(errstr)

    def get_results(self, resultspath):
        '''Moves work from boinc results directory to results path.'''
        # TODO: Think about writing a boinc assimilator that does this stuff.
        all_boinc_results = os.listdir(self.boinc_results_path)
        all_boinc_results = [ f for f in all_boinc_results if '_' in f ]
        my_boinc_results = [ f for f in all_boinc_results
                                if int(f.split('_')[0]) == self.uniqueid and
                                'no_output_files' not in f]

        #build up a set of all the state_wunumbers 
        jobfiles = {}
        for fname in my_boinc_results:
            fname_parts = fname.split('_')    
            fname_path = os.path.join(self.boinc_results_path,fname)
            # fnames look like uniqueid_statenumber_wunumber_filenumber
            # split on underscores  0        1          2        3
            # thus key is statenumber_wunumber
            key = fname_parts[1]+'_'+fname_parts[2]
            if key in jobfiles:
                jobfiles[key].append(fname_path)
            else:
                jobfiles[key] = [fname_path]

        # XXX: This is completely dependent on the result template file.
        #      We should be parsing the result template file for this stuff.                
        ending_to_filename = [ 'results.dat', 'reactant.con', 'saddle.con', 
                'product.con', 'mode.dat' ]

        keys = jobfiles.keys()

        for key in keys:
            filepaths = jobfiles[key]
            num_files = len(filepaths)
            if num_files != 5:
                logger.warning("got %i out of 5 files for %s" % (num_files, key))
                jobfiles.pop(key)
                continue

            try:
                os.makedirs(os.path.join(resultspath, key))
            except OSError, (errno, strerrno):
                if errno == 31:
                    logger.warning(strerrno)
                    break
                else:
                    raise

            for filepath in filepaths:
                path, filename = os.path.split(filepath)
                ending = int(filename[-1])
                filename = ending_to_filename[ending]
                shutil.move(filepath, os.path.join(resultspath, key, filename))

        return self.unbundle(resultspath)

class Local(Communicator):
    def __init__(self, scratchpath, client, ncpus, bundle_size):
        Communicator.__init__(self, scratchpath, bundle_size)

        #number of cpus to use
        self.ncpus = ncpus
        #path to the saddle search executable
        self.client = client

        self.searchlist = []

        import atexit
        #don't let clients hang around if the script dies
        atexit.register(self.cleanup)
        
    def cleanup(self):
        '''Kills the running eonclients.'''
        import signal
        for search in self.searchlist:
            p = search[0]
            try:
                os.kill(p.pid, signal.SIGKILL)
            except:
                pass

    def get_results(self, resultspath):
        '''Moves work from scratchpath to results path.'''
        jobdirs = [ d for d in os.listdir(self.scratchpath) 
                    if os.path.isdir(os.path.join(self.scratchpath,d)) ]

        for jobdir in jobdirs:
            dest_dir = os.path.join(resultspath, jobdir)
            shutil.move(os.path.join(self.scratchpath,jobdir), dest_dir)
        return self.unbundle(resultspath)

        #jobdirs = [ os.path.join(self.scratchpath, d) for d in os.listdir(self.scratchpath) 
        #                if os.path.isdir(os.path.join(self.scratchpath,d)) ]
        #results = []
        #for jobdir in self.unbundle(jobdirs):
        #    shutil.move(jobdir, os.path.join(resultspath, os.path.split(jobdir)[1]))
        #    results.append(os.path.split(jobdir)[1])
        #return results

    def check_search(self, search):
        p, jobpath = search
        if p.returncode == 0:
            logger.info('saddle search finished in %s' % jobpath)
            return True
        else:
            stdout, stderr = p.communicate()
            errmsg = "saddle search failed in %s: %s" % (jobpath, stderr)
            logger.warning(errmsg)

    def submit_searches(self, jobpaths):
        '''Run up to ncpu number of clients to process the work in jobpaths.
           The job directories are moved to the scratch path before the calculcation
           is run. This method doesn't return anything.'''
        for jobpath in self.make_bundles(jobpaths):
            #move the job directory to the scratch directory
            jobdir = os.path.split(jobpath)[1]
            destpath = os.path.join(self.scratchpath, jobdir)
            shutil.move(jobpath, destpath)
            #update jobpath to be in the scratch directory
            jobpath = destpath
            fstdout = open(os.path.join(jobpath, "stdout.dat"),'w')
            p = subprocess.Popen(self.client,cwd=jobpath,
                    stdout=fstdout, stderr=subprocess.PIPE)
            commands.getoutput("renice -n 20 -p %d" % p.pid)
            self.searchlist.append((p,jobpath))

            while len(self.searchlist) == self.ncpus:
                for i in range(len(self.searchlist)):
                    p = self.searchlist[i][0]
                    retval = p.poll()
                    if retval == None:
                        continue
                    else:
                        self.check_search(self.searchlist[i])
                        self.searchlist.pop(i)
                        break
                sleep(0.1)

        #wait for everything to finish
        for search in self.searchlist:
            p = search[0]
            p.wait()
            self.check_search(search)

    def cancel_state(self, state):
        return 0

    def get_queue_size(self):
        return 0


class ARC(Communicator):

    def __init__(self, scratchpath, bundle_size=1):
        Communicator.__init__(self, scratchpath, bundle_size)

        try:
            import arclib
            self.arc = arclib
        except ImportError:
            raise CommunicatorError("ARCLib can't be imported. Check if PYTHONPATH is set correctly")

        self.jobsfilename = os.path.join(self.scratchpath, "jobs.txt")


    def create_wrapper_script(self):
        '''Create a wrapper script to execute a job.  Return path to script.'''
        
        s = """
        #!/bin/sh

        tar jxvf $1.tar.bz2
        cd $1
        Client
        cd $HOME
        tar jcvf $1.tar.bz2  $1
        """
        script_path = os.path.join(self.scratchpath, 'wrapper.sh')
        try:
            f = open(script_path, "w")
            f.write(s)
            f.close()
        except Exception as msg:
            raise CommunicatorError("Can't create wrapper script: %s" % msg)

        return script_path


    def create_tarball(self, src, dest):
        """Pack directory 'src' into tar.bz2 file 'dest', makeing sure it will unopack into
           a dir called basename(src), rather than path/to/src"""

        # Remove trailing '/'; it would cause trouble with
        # os.path.dirname() and os.path.basename()
        if src[-1] == '/':
            src = src[:-1]

        dirname = os.path.dirname(src)
        basename = os.path.basename(src)

        cwd = os.getcwd()
        try:
            os.chdir(dirname)
            tarball = tarfile.open(dest, 'w:bz2')
            tarball.add(basename)
            tarball.close()
        finally:
            os.chdir(cwd)


    def create_job(self, job_path, wrapper_path):
        '''Prepare a job who's inputfiles are found in 'job_path'.
           Return pre-processed xRSL code and job name.'''

        # Remove trailing '/'; it would cause trouble with os.path.basename()
        if job_path[-1] == '/':
            job_path = job_path[:-1]

        basename = os.path.basename(job_path)
        tarball_path = os.path.join(self.scratchpath, basename + ".tar.bz2")
        self.create_tarball(job_path, tarball_path)

        s = "&"
        s += "(executable=./%s)" % os.path.basename(wrapper_path)
        s += "(arguments=%s)" % basename

        s += "(inputFiles="
        s += "(%s %s)" % (os.path.basename(wrapper_path), wrapper_path)
        s += "(%s %s)" % (os.path.basename(tarball_path), tarball_path)
        s += ")"

        s += "(outputFiles="
        s += "(%s '')" % os.path.basename(tarball_path)
        s += ")"

        s += "(stdout=stdout)"
        s += "(stderr=stderr)"

        s += "(runTimeEnvironment=APPS/CHEM/EON2)"

        jobname = "eon-%s" % basename
        s += "(jobName=%s)" % jobname

        logger.debug("xrsl: " + s)

        return self.arc.Xrsl(s), jobname


    def submit_searches(self, jobpaths):
        '''Throws CommunicatorError if fails.'''

        try:
            c = self.arc.Certificate(self.arc.PROXY)
        except self.arc.CertificateError as msg:
            raise CommunicatorError(msg)

        if c.IsExpired():
            raise CommunicatorError("Grid proxy has expired!")

        logger.info("Grid proxy is valid for " + c.ValidFor())

        wrapper_path = self.create_wrapper_script()

        qi = self.arc.GetQueueInfo()
        f = open(self.jobsfilename, "w")

        for job_path in self.make_bundles(jobpaths):
            xrsl, jobname = self.create_job(job_path, wrapper_path)
            targets = self.arc.ConstructTargets(qi, xrsl)
            targetsleft = self.arc.PerformStandardBrokering(targets)
            jobid = self.arc.SubmitJob(xrsl, targetsleft)

            self.arc.AddJobID(jobid, jobname) # Why is this needed?
            f.write(jobid + '\n')
            logger.debug("jobId: " + jobid)
        f.close()


    def get_results(self, resultspath):
        '''Returns a list of directories containing the results.'''
        return []

    def get_queue_size(self):
        '''Returns the number of items waiting to run in the queue.'''
        return 0

    def cancel_state(self, statenumber):
        '''Returns the number of workunits that were canceled.'''
        print 'statenumber =', statenumber
        raise NotImplementedError()



if __name__=='__main__':
    import sys
    if len(sys.argv) < 4:
        print '%s: scratchdir jobdir resultsdir' % sys.argv[0]
        sys.exit(1)
    c = Local('eonclient', 2, sys.argv[1]) 
    jobdirs = [ os.path.join(sys.argv[2], d) for d in os.listdir(sys.argv[2]) ]

    c.submit_searches(jobdirs)
    c.get_results(sys.argv[3])
