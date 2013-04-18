### subroutines for creation of subdirectories & clean-up

import os

def mklocaltmp(odir, site):
    if site.batch:
        s = site.submitdir
        job = site.jobid
    else:
        s = os.getcwd()
        job = ''
    if odir is None or len(odir)==0:
        p = os.popen('mktemp -d '+s+'/qe'+job+'_XXXXX', 'r')
        tdir = p.readline().strip()
        p.close()
    else:
        if odir[0]=='/':
            tdir = odir
        else:
            tdir = s+'/'+odir
        os.system('mkdir -p '+tdir)
    return tdir

def mkscratch(localtmp, site):
    if site.batch:
        pernodeexec = site.perHostMpiExec
        job = site.jobid
    else:
        pernodeexec = ''
        job = ''
    p = os.popen('mktemp -d '+site.scratch+'/qe'+job+'_XXXXX', 'r')
    tdir = p.readline().strip()
    p.close()
    if pernodeexec!='':
        cdir = os.getcwd()
        os.chdir(localtmp)
        os.system(pernodeexec + ' mkdir -p '+tdir)
        os.chdir(cdir)
    return tdir

def cleanup(tmp, scratch, removewf, calc, site):
    try:
        calc.stop()
    except:
        pass
    if site.batch:
        pernodeexec = site.perHostMpiExec
    else:
        pernodeexec = ''
    if removewf:
        os.system('rm -r '+scratch+'/*wfc* 2>/dev/null')
    os.system('cp -r '+scratch+' '+tmp)
    cdir = os.getcwd()
    os.chdir(tmp)
    os.system(pernodeexec + ' rm -r '+scratch+' 2>/dev/null')
    os.chdir(cdir)

def getsubmitorcurrentdir(site):
    s = site.submitdir
    if s!='':
        return s
    else:
        return os.getcwd()
