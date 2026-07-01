
'''
Con(figuration) i/o library

CON I/O is implemented on top of the modern **readcon** Python bindings
(readcon-core): vectorized ``coords_array`` / velocities / forces, frame
metadata, first-frame and streaming iterators, and optional gzip writes.
Callers still use :func:`loadcon` / :func:`savecon` and an :class:`eon.atoms.Atoms`
view; optional fields on ``Atoms`` (``v``, ``f``, ``energy_per_atom``,
``frame_metadata``, ``fixed_axes``) carry full-frame payload when present.
'''
import configparser
#from io import BytesIO as StringIO
from io import StringIO
import logging
logger = logging.getLogger('io')
import numpy
import os

import pickle as pickle
import readcon

from eon import atoms
from eon.config import config

def save_prng_state():
    state = numpy.random.get_state()
    fh = open('prng.pkl', 'wb')
    pickle.dump(state, fh, pickle.HIGHEST_PROTOCOL)

def get_prng_state():
    fh = open('prng.pkl', 'rb')
    state = pickle.load(fh)
    numpy.random.set_state(state)

def length_angle_to_box(boxlengths, angles):
    box = numpy.zeros( (3,3) )
    angles *= numpy.pi/180.0
    box[0][0] = 1.0
    box[1][0] = numpy.cos(angles[0])
    box[1][1] = numpy.sin(angles[0])
    box[2][0] = numpy.cos(angles[1])
    box[2][1] = (numpy.cos(angles[2])-box[1][0]*box[2][0])/box[1][1]
    box[2][2] = numpy.sqrt(1.0-box[2][0]**2-box[2][1]**2)
    box[0,:]*=boxlengths[0]
    box[1,:]*=boxlengths[1]
    box[2,:]*=boxlengths[2]
    return box

def box_to_length_angle(box):
    lengths = numpy.zeros(3)
    lengths[0] = numpy.linalg.norm(box[0,:])
    lengths[1] = numpy.linalg.norm(box[1,:])
    lengths[2] = numpy.linalg.norm(box[2,:])
    angles = numpy.zeros(3)
    angles[0] = numpy.arccos(numpy.dot(box[0,:]/lengths[0],box[1,:]/lengths[1]))
    angles[1] = numpy.arccos(numpy.dot(box[0,:]/lengths[0],box[2,:]/lengths[2]))
    angles[2] = numpy.arccos(numpy.dot(box[1,:]/lengths[1],box[2,:]/lengths[2]))
    angles *= 180.0/numpy.pi
    return lengths, angles


def _frame_metadata_dict(frame):
    """Return a plain dict copy of frame.metadata when present."""
    meta = getattr(frame, 'metadata', None)
    if meta is None:
        return None
    try:
        return dict(meta)
    except TypeError:
        return None


def _frame_to_atoms(frame):
    """Convert a :class:`readcon.ConFrame` to an eon :class:`atoms.Atoms`.

    Uses vectorized ``coords_array`` / ``velocities_array`` / ``forces_array`` /
    ``energies_array`` when available (readcon >= 0.13), falling back to
    per-atom Python attributes for older wheels.
    """
    num_atoms = len(frame)
    a = atoms.Atoms(num_atoms)
    boxlengths = numpy.asarray(frame.cell, dtype=float)
    boxangles = numpy.asarray(frame.angles, dtype=float)
    a.box = length_angle_to_box(boxlengths, boxangles)

    # Positions: prefer contiguous ndarray from the bindings.
    coords = None
    if hasattr(frame, 'coords_array'):
        try:
            coords = numpy.asarray(frame.coords_array(), dtype=float)
        except Exception:
            coords = None
    frame_atoms = list(frame.atoms)
    if coords is not None and coords.shape == (num_atoms, 3):
        a.r[:] = coords
    else:
        for i, atom in enumerate(frame_atoms):
            a.r[i] = [atom.x, atom.y, atom.z]

    fixed_axes = numpy.zeros((num_atoms, 3), dtype=bool)
    for i, atom in enumerate(frame_atoms):
        a.names[i] = atom.symbol
        mass = atom.mass
        a.mass[i] = float(mass) if mass is not None else 0.0
        fixed = list(atom.fixed) if atom.fixed is not None else [False, False, False]
        if len(fixed) != 3:
            fixed = [bool(fixed[0])] * 3 if fixed else [False, False, False]
        fixed_axes[i] = fixed
        # Legacy free flag: free if no axis is fixed.
        a.free[i] = 0 if any(fixed) else 1
    a.fixed_axes = fixed_axes

    # Optional velocity / force / energy blocks.
    if hasattr(frame, 'velocities_array'):
        try:
            vel = frame.velocities_array()
            if vel is not None:
                a.v = numpy.asarray(vel, dtype=float)
        except Exception:
            pass
    if a.v is None and frame_atoms and getattr(frame_atoms[0], 'has_velocity', False):
        a.v = numpy.zeros((num_atoms, 3))
        for i, atom in enumerate(frame_atoms):
            a.v[i] = [
                atom.vx if atom.vx is not None else 0.0,
                atom.vy if atom.vy is not None else 0.0,
                atom.vz if atom.vz is not None else 0.0,
            ]

    if hasattr(frame, 'forces_array'):
        try:
            frc = frame.forces_array()
            if frc is not None:
                a.f = numpy.asarray(frc, dtype=float)
        except Exception:
            pass
    if a.f is None and frame_atoms and getattr(frame_atoms[0], 'has_forces', False):
        a.f = numpy.zeros((num_atoms, 3))
        for i, atom in enumerate(frame_atoms):
            a.f[i] = [
                atom.fx if atom.fx is not None else 0.0,
                atom.fy if atom.fy is not None else 0.0,
                atom.fz if atom.fz is not None else 0.0,
            ]

    if hasattr(frame, 'energies_array'):
        try:
            eng = frame.energies_array()
            if eng is not None:
                a.energy_per_atom = numpy.asarray(eng, dtype=float)
        except Exception:
            pass
    if a.energy_per_atom is None and frame_atoms and getattr(
            frame_atoms[0], 'has_energy', False):
        a.energy_per_atom = numpy.array(
            [atom.energy if atom.energy is not None else 0.0 for atom in frame_atoms],
            dtype=float,
        )

    a.frame_metadata = _frame_metadata_dict(frame)
    return a


def _read_all_frames_path(path):
    """All frames from a path (``read_all_frames`` or ``read_con``)."""
    reader = getattr(readcon, 'read_all_frames', None) or readcon.read_con
    return reader(path)


def loadcons(filename):
    """Load every frame in a multi-frame ``.con`` / ``.convel`` / ``.con.gz``."""
    return [_frame_to_atoms(f) for f in _read_all_frames_path(filename)]


def iter_cons(filename):
    """Yield :class:`atoms.Atoms` for each frame (streaming via ``iter_con``)."""
    if hasattr(readcon, 'iter_con'):
        for frame in readcon.iter_con(filename):
            yield _frame_to_atoms(frame)
    else:
        for frame in _read_all_frames_path(filename):
            yield _frame_to_atoms(frame)


def loadposcars(filename):
    filein = open(filename, 'r')
    p = []
    while True:
        try:
            p.append(loadposcar(filein))
        except:
            return p


def loadcon(filein, reset=True):
    """
    Load the **first** frame of a con file.

    ``filein`` may be a filename or a file-like object. Paths use
    :func:`readcon.read_first_frame` (no full multi-frame parse). File-like
    objects use :func:`readcon.read_con_string`.
    """
    if hasattr(filein, 'readline'):
        content = filein.read()
        if reset:
            try:
                filein.seek(0)
            except Exception:
                pass
        frames = readcon.read_con_string(content)
        if not frames:
            raise IOError("No frames found in con data")
        return _frame_to_atoms(frames[0])

    path = filein
    if hasattr(readcon, 'read_first_frame'):
        try:
            return _frame_to_atoms(readcon.read_first_frame(path))
        except Exception as exc:
            # Fall back to full read for odd extensions / older edge cases.
            logger.debug("read_first_frame failed for %s (%s); using read_con",
                         path, exc)
    frames = readcon.read_con(path)
    if not frames:
        raise IOError("No frames found in con data")
    return _frame_to_atoms(frames[0])


def loadcon_as_ase(filename):
    """Load CON path as a list of ASE Atoms via readcon (requires ``ase``)."""
    if hasattr(readcon, 'read_con_as_ase'):
        return readcon.read_con_as_ase(filename)
    raise RuntimeError("readcon.read_con_as_ase is not available in this wheel")


def _atoms_to_frame(p, metadata=None):
    """Convert an eon :class:`atoms.Atoms` to a :class:`readcon.ConFrame`.

    Propagates optional velocities, forces, per-atom energies, per-axis fixed
    flags, and frame metadata when present on ``p``.
    """
    lengths, angles = box_to_length_angle(p.box)
    atom_list = []
    n = len(p)
    use_fixed_axes = (
        getattr(p, 'fixed_axes', None) is not None
        and numpy.asarray(p.fixed_axes).shape == (n, 3)
    )
    use_v = getattr(p, 'v', None) is not None and numpy.asarray(p.v).shape == (n, 3)
    use_f = getattr(p, 'f', None) is not None and numpy.asarray(p.f).shape == (n, 3)
    use_e = (
        getattr(p, 'energy_per_atom', None) is not None
        and numpy.asarray(p.energy_per_atom).shape == (n,)
    )
    for i in range(n):
        if use_fixed_axes:
            fixed = [bool(p.fixed_axes[i, 0]), bool(p.fixed_axes[i, 1]),
                     bool(p.fixed_axes[i, 2])]
        else:
            fixed = [p.free[i] == 0] * 3
        kwargs = dict(
            symbol=p.names[i],
            x=float(p.r[i][0]),
            y=float(p.r[i][1]),
            z=float(p.r[i][2]),
            fixed=fixed,
            atom_id=i + 1,
            mass=float(p.mass[i]),
        )
        if use_v:
            kwargs['vx'] = float(p.v[i][0])
            kwargs['vy'] = float(p.v[i][1])
            kwargs['vz'] = float(p.v[i][2])
        if use_f:
            kwargs['fx'] = float(p.f[i][0])
            kwargs['fy'] = float(p.f[i][1])
            kwargs['fz'] = float(p.f[i][2])
        if use_e:
            kwargs['energy'] = float(p.energy_per_atom[i])
        atom_list.append(readcon.Atom(**kwargs))

    meta = metadata
    if meta is None and getattr(p, 'frame_metadata', None) is not None:
        meta = dict(p.frame_metadata)
    frame_kwargs = dict(
        cell=[float(x) for x in lengths],
        angles=[float(x) for x in angles],
        atoms=atom_list,
        prebox_header=["Generated by eOn", ""],
    )
    # metadata= supported on modern ConFrame; omit if constructor rejects it.
    if meta:
        try:
            return readcon.ConFrame(metadata=meta, **frame_kwargs)
        except TypeError:
            pass
    return readcon.ConFrame(**frame_kwargs)


def _write_con_path(path, frames, precision=6, compression=None, canonical=None):
    """Call :func:`readcon.write_con` with optional kwargs only if supported."""
    kwargs = {'precision': precision}
    if compression is not None:
        kwargs['compression'] = compression
    # canonical= added in newer wheels; ignore if unsupported.
    if canonical is not None:
        try:
            readcon.write_con(path, frames, canonical=canonical, **kwargs)
            return
        except TypeError:
            pass
    readcon.write_con(path, frames, **kwargs)


def savecon(fileout, p, w='w', precision=6, compression=None, canonical=None,
            metadata=None):
    """
    Save a con file using readcon.

    Parameters
    ----------
    fileout : str or file-like
        Path or writable object.
    p : atoms.Atoms
        Configuration to write.
    w : {'w', 'a'}
        Write or append (append rewrites multi-frame via readcon for paths).
    precision : int
        Decimal precision for coordinates (passed to readcon).
    compression : {None, 'gzip', 'none'}
        Force compression mode for path writes (default: auto from ``.gz``).
    canonical : bool or None
        If True and supported by the installed readcon, write canonical CON.
    metadata : dict or None
        Optional frame metadata override (else ``p.frame_metadata``).
    """
    frame = _atoms_to_frame(p, metadata=metadata)
    if hasattr(fileout, 'write'):
        text_kwargs = {'precision': precision}
        try:
            if canonical is not None:
                text = readcon.write_con_string([frame], canonical=canonical,
                                                **text_kwargs)
            else:
                text = readcon.write_con_string([frame], **text_kwargs)
        except TypeError:
            text = readcon.write_con_string([frame], precision=precision)
        fileout.write(text)
        return

    path = fileout
    if w == 'a' and os.path.exists(path):
        existing = list(_read_all_frames_path(path))
        existing.append(frame)
        _write_con_path(path, existing, precision=precision,
                        compression=compression, canonical=canonical)
    else:
        _write_con_path(path, [frame], precision=precision,
                        compression=compression, canonical=canonical)


def savecons(fileout, atoms_list, precision=6, compression=None, canonical=None):
    """Write a multi-frame CON trajectory in one readcon call."""
    frames = [_atoms_to_frame(p) for p in atoms_list]
    if hasattr(fileout, 'write'):
        try:
            text = readcon.write_con_string(frames, precision=precision)
        except TypeError:
            text = readcon.write_con_string(frames)
        fileout.write(text)
        return
    _write_con_path(fileout, frames, precision=precision,
                    compression=compression, canonical=canonical)


def load_mode(modefilein):
    '''
    Reads a mode.dat file into an N by 3 numpy array
        modefilein: may be either a file-like object of a filename
    '''
    if hasattr(modefilein, 'readline'):
        f = modefilein
    else:
        f = open(modefilein, 'r')
    if len(f.readline().split()) == 3:
        f.seek(0);
    lines = f.readlines()
    mode = []
    for line in lines:
        l = line.strip().split()
        for j in range(3):
            mode.append(float(l[j]))
    return numpy.array(mode).reshape(len(mode)//3, 3)

def save_mode(modefileout, displace_vector):
    '''
    Saves an Nx3 numpy array into a mode.dat file.
        modefileout:     may be either a filename or file-like object
        displace_vector: the mode (Nx3 numpy array)
    '''
    if hasattr(modefileout, 'write'):
        f = modefileout
    else:
        f = open(modefileout, 'w')
    for i in range(len(displace_vector)):
        f.write("%.3f %.3f %.3f\n" % (displace_vector[i][0],
            displace_vector[i][1], displace_vector[i][2]))


def save_results_dat(fileout, results):
    '''
    Saves a results.dat file from a dictionary
    '''
    if hasattr(fileout, 'write'):
        f = fileout
    else:
        f = open(fileout, 'w')

    for key in results:
        #print(results[key], key, con)  #GH: this made no sense to me - replaced with the following
        f.write(results[key], key)

def modify_config(config_path, changes):
    parser = configparser.ConfigParser()
    parser.read(config.config_path)
    for change in changes:
        parser.set(*change)
    config_str_io = StringIO()
    parser.write(config_str_io)
    config_str_io.seek(0)
    return config_str_io

def parse_results(filein):
    '''
    Reads a results.dat file and gives a dictionary of the values contained therein
    '''
    if hasattr(filein, 'readline'):
        f = filein
        f.seek(0)
    else:
        f = open(filein)
    results = {}
    for line in f:
        line = line.split()
        if len(line) < 2:
            continue
        if '.' in line[0]:
            try:
                results[line[1]] = float(line[0])
            except ValueError:
                logger.warning("Couldn't parse float in results.dat: %s", line)
        else:
            try:
                results[line[1]] = int(line[0])
            except ValueError:
                try:
                    results[line[1]] = line[0]
                except ValueError:
                    logger.warning("Couldn't parse string in results.dat: %s", line)

    return results

def loadposcar(filein):
    '''
    Load the POSCAR file named filename and returns an atoms object
    '''
    if hasattr(filein, 'readline'):
        f = filein
    else:
        f = open(filein, 'r')
    # Line 1: Atom types
    AtomTypes = f.readline().split()
    # Line 2: scaling of coordinates
    scale = float(f.readline())
    # Lines 3-5: the box
    box = numpy.zeros((3, 3))
    for i in range(3):
        line = f.readline().split()
        box[i] = numpy.array([float(line[0]), float(line[1]), float(line[2])]) * scale
    # Line 6: number of atoms of each type.
    line = f.readline().split()
    NumAtomsPerType = []
    for l in line:
        NumAtomsPerType.append(int(l))
    # Now have enough info to make the atoms object.
    num_atoms = sum(NumAtomsPerType)
    p = atoms.Atoms(num_atoms)
    # Fill in the box.
    p.box = box
    # Line 7: selective or cartesian
    sel = f.readline()[0]
    selective_flag = (sel == 's' or sel == 'S')
    if not selective_flag:
        car = sel
    else:
        car = f.readline()[0]
    direct_flag = not (car == 'c' or car == 'C' or car == 'k' or car == 'K')
    atom_index = 0
    for i in range(len(NumAtomsPerType)):
        for j in range(NumAtomsPerType[i]):
            p.names[atom_index] = AtomTypes[i]
            line = f.readline().split()
            if(selective_flag):
                assert len(line) >= 6
            else:
                assert len(line) >= 3
            pos = line[0:3]
            if selective_flag:
                sel = line[3:7]
                if sel[0] == 'T' or sel[0] == 't':
                    p.free[atom_index] = 1
                elif sel[0] == 'F' or sel[0] == 'f':
                    p.free[atom_index] = 0
            p.r[atom_index] = numpy.array([float(q) for q in pos])
            if direct_flag:
                p.r[atom_index] = numpy.dot(p.r[atom_index], p.box)
            else:
                p.r[atom_index] *= scale
            atom_index += 1
    return p


def saveposcar(fileout, p, w='w', direct = False):
    '''
    Save a POSCAR
        fileout: name to save it under
        point:    atoms object to save
        w:        write/append flag
    '''
    if hasattr(fileout, 'write'):
        poscar = fileout
    else:
        poscar = open(fileout, w)
    atom_types = []
    num_each_type = {}
    for name in p.names:
        if not name in atom_types:
            atom_types.append(name)
            num_each_type[name] = 1
        else:
            num_each_type[name] += 1
    poscar.write(" ".join(atom_types)+'\n') #Line 1: Atom type
    poscar.write("1.0\n") #Line 2: scaling
    for i in range(3):
        poscar.write(" ".join(['%20.14f' % s for s in p.box[i]])+'\n')  #lines 3-5: box
    poscar.write(" ".join(atom_types)+'\n') #Line 6: Atom type
    poscar.write(" ".join(['%s' % num_each_type[key] for key in atom_types])+'\n')
    poscar.write('Selective Dynamics\n') #line 7: selective dynamics
    if direct:
        poscar.write('Direct\n')  #line 8 cartesian coordinates
        ibox = numpy.linalg.inv(numpy.array(p.box))
        p.r = numpy.dot(p.r, ibox)
    else:
        poscar.write('Cartesian\n') #line 8 cartesian coordinates
    for i in range(len(p)):
            posline = " ".join(['%20.14f' % s for s in p.r[i]]) + " "
            for j in range(3):
                if(p.free[i]):
                    posline+='   T'
                else:
                    posline+='   F'
            poscar.write(posline+'\n')


from configparser import ConfigParser as SCP
class ini(SCP):

    def __init__(self, filenames):
        self.loaded = False
        self.filenames = filenames
        SCP.__init__(self)

    def read(self):
        self.loaded = True
        SCP.read(self, self.filenames)

    def get(self, section, option, default="ini_no_default", **kwargs):
#    def get(self, section, option, default="ini_no_default"):
        if not self.loaded:
            self.read()
        try:
            SCP.read(self, self.filenames)
            #value = SCP.get(self, section, option, raw=True, **kwargs)
            value = SCP.get(self, section, option, **kwargs)
            #value = SCP.get(self, section, option, raw=True)
            #value = SCP.get(self, section, option)
        except:
            if default == "ini_no_default":
                raise NameError("Section or option missing, no default specified")
            return default
        try:
            return int(value)
        except ValueError:
            pass
        try:
            return float(value)
        except ValueError:
            pass
        if value.lower() == 'true':
            return True
        if value.lower() == 'false':
            return False
        return value

    def getint(self, *args):
        raise NotImplementedError("Use the get function with this ConfigParser wrapper")
    def getfloat(self, *args):
        raise NotImplementedError("Use the get function with this ConfigParser wrapper")
    def getboolean(self, *args):
        raise NotImplementedError("Use the get function with this ConfigParser wrapper")

    def set(self, section, option, value):
        if section not in self.sections():
            self.add_section(section)
        SCP.set(self, section, option, str(value))
        if type(self.filenames) == str:
            name = self.filenames
        else:
            name = self.filenames[-1]
#        configfile = open(name, 'wb')
        configfile = open(name, 'w')
        self.write(configfile)
        configfile.close()


class Dynamics:
    """ The Dynamics class handles I/O for the dynamics.txt file of an aKMC simulation. """

    def __init__(self, filename):
        self.filename = filename
        if not os.path.exists(filename):
            f = open(self.filename, 'w')
            header = "%12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s\n" % ('step-number', 'reactant-id', 'process-id', 'product-id', 'step-time', 'total-time', 'barrier', 'rate', 'energy')
            f.write(header)
            f.write("-" * len(header))
            f.write("\n")
            f.close()
            self.next_step = 0

        # read last lines of the file to determine iteration nr
        else:
            f = open(self.filename,'r')
            f.seek(0,2)	#seek to EOF
            fsize = f.tell()
            # seek 1024 bytes back (or to beginning of file if fsize < 1024 )
            # last line must be contained in this block
            f.seek( max( fsize - 1024 , 0 ) , 0)
            lines = f.readlines()
            self.next_step = int ( lines[-1].split()[0] ) + 1 # determine iteration nr of next step

    def append(self, reactant_id, process_id, product_id, step_time, total_time, barrier, rate, energy):
        f = open(self.filename, 'a')
        f.write("%12d  %12d  %12d  %12d  %12e  %12e  %12f  %12e  %12f\n" % (self.next_step, reactant_id, process_id, product_id, step_time, total_time, barrier, rate, energy))
        f.close()
        self.next_step += 1

    def append_sb(self, reactant_id, process_id, product_id, step_time, total_time, basin_id, rate, energy):
        f = open(self.filename, 'a')
        f.write("%12d  %12d  %12d  %12d  %12e  %12e  %12d  %12e  %12f\n" % (self.next_step, reactant_id, process_id, product_id, step_time, total_time, basin_id, rate, energy))
        f.close()
        self.next_step += 1

    def get(self):
        f = open(self.filename, 'r')
        lines = f.readlines()[2:]
        f.close()
        data = []
        for line in lines:
            split = line.split()
            data.append({"reactant":    int(split[1]),
                         "process":     int(split[2]),
                         "product":     int(split[3]),
                         "steptime":    float(split[4]),
                         "totaltime":   float(split[5]),
                         "barrier":     float(split[6]),
                         "prefactor":   float(split[7])})
        return data

def load_potfiles(pot_dir):
    ret = {}
    if os.path.isdir(pot_dir):
        for name in os.listdir(pot_dir):
            if os.path.isdir(name):
                continue
            a = open(os.path.join(pot_dir, name), 'r')
            b = StringIO("".join(a.readlines()))
            c = os.stat(os.path.join(pot_dir, name)).st_mode
            ret[name] = (b,c)
    return ret

class TableException(Exception):
    pass

class Table:
    """
    A class that provides a nice io abstraction for table like data.  The data
    is saved in a pretty printed format. Also provides nice data retrival
    methods.

    >>> t = Table("sample.tbl", ['id', 'name', 'age' ])
    >>> t.eagerwrite = False
    >>> t.add_row({'id':0,'name':"Sam","age":24})
    >>> t.add_row({'id':1,'name':"David","age":50})
    >>> t.add_row({'id':2,'name':"Anna","age":21})
    >>> t #doctest: +NORMALIZE_WHITESPACE
        id name  age
        -- ----- ---
        0  Sam   24
        1  David 50
        2  Anna  21

    Rows can be accessed directly:
    >>> t.rows[1] #doctest: +SKIP
        {'age': 50, 'id': 1, 'name': 'David'}
    >>> t.max_value('age') #doctest: +NORMALIZE_WHITESPACE
        50
    >>> t.min_row('age') #doctest: +NORMALIZE_WHITESPACE +SKIP
        {'age': 21, 'id': 2, 'name': 'Anna'}
    >>> sorted(t.min_row('id').items()) #doctest: +NORMALIZE_WHITESPACE
        [('age', 24), ('id', 0), ('name', 'Sam')]
    >>> len(t) #doctest: +NORMALIZE_WHITESPACE
        3
    >>> sum(t.getcolumn('age')) #doctest: +NORMALIZE_WHITESPACE
        95
    >>> t.write() #doctest: +SKIP

    The table can be loaded from disk without specifying columns. This is
    slightly unsafe because the columns can't be checked, but it could cut down
    on the verbosity in some places.
    >>> t2 = Table("sample.tbl") #doctest: +SKIP
"""

    #XXX: This is the number of digits that a floating point number gets
    #     serialized with. Should it be some sort of config option?
    #     Or is there just a good default?

    def __init__(self, filename, columns=None, overwrite=False):
        self.filename = filename
        self.columns = columns
        self.rows = []
        self.columntypes = {}
        self.columnwidths = {}
        self.initialized = False
        self.overwrite = overwrite

        self.floatprecision = 6
        self.eagerwrite = True

    def init(self):
        """Checks to see if self.filename exists. If it does self.rows
           will be initialized from disk."""
        self.initialized = True
        if os.path.isfile(self.filename) and not self.overwrite:
            self.read(self.filename)
        else:
            if self.columns is None:
                raise TableException("Columns are not optional for new tables")

            for c in self.columns:
                self.columnwidths[c] = len(c)

    def read(self, filename):
        self.eagerwrite = False
        f = open(self.filename, "r")
        filecolumns = f.readline().split()
        if self.columns != None:
            if filecolumns != self.columns:
                raise TableException("Column name mismatch: %s" % filename)
        else:
            self.columns = filecolumns

        for c in self.columns:
            self.columnwidths[c] = len(c)

        # skip comment line
        f.readline()

        for line in f:
            fields = line.split()
            row = {}
            coli = 0
            for field in fields:
                try:
                    field = int(field)
                except ValueError:
                    try:
                        field = float(field)
                    except ValueError:
                        field = field.strip()
                        pass
                row[self.columns[coli]] = field
                coli += 1
            self.add_row(row)
        f.close()
        self.eagerwrite = True

    def __repr__(self):
        if not self.initialized:
            self.init()
        f = StringIO()
        self.writefilehandle(f)
        return f.getvalue()

    def __len__(self):
        if not self.initialized:
            self.init()
        return len(self.rows)

    def __iter__(self):
        for row in self.rows:
            yield row

    def write(self):
        if not self.initialized:
            self.init()
        f = open(self.filename, "w")
        #print("into table write: ",self.filename)
        self.writefilehandle(f)
        f.close()

    def writefilehandle(self, filehandle):
        f = filehandle
        line = ' '.join([ "%-*s"%(self.columnwidths[c], c) for c in self.columns ])
        f.write(line+"\n")

        line = ''
        for c in self.columns:
            line += '-'*self.columnwidths[c]+' '
        f.write(line+'\n')

        for row in self.rows:
            line = ""
            for c in self.columns:
                if self.columntypes[c] == float:
                    line += "%#-*.*G " % (self.columnwidths[c],self.floatprecision,
                                         row[c])
                else:
                    line += "%-*s " % (self.columnwidths[c],row[c])
            f.write(line+"\n")

    def add_row(self, row):
        if not self.initialized:
            self.init()
        mismatched_columns = set(self.columns).symmetric_difference(set(row.keys()))
        if len(mismatched_columns) != 0:
            raise TableException("Mismatched columns %s" % str(mismatched_columns))

        if len(self.rows) == 0:
            for c in row:
                self.columntypes[c] = type(row[c])
        else:
            for c in row:
                if type(row[c]) != self.columntypes[c]:
                    raise TableException("Type mismatch for column %s" % c)

        for c in row:
            if self.columntypes[c] == float:
                self.columnwidths[c] = max(self.columnwidths[c], self.floatprecision+5)
            else:
                self.columnwidths[c] = max(self.columnwidths[c], len(str(row[c])))

        self.rows.append(row)
        if self.eagerwrite:
            self.write()

    def delete_row(self, column, value):
        if not self.initialized:
            self.init()
        rows_to_delete = []
        for row in self.rows:
            if row[column] == value:
                rows_to_delete.append(row)
#        map(self.rows.remove, rows_to_delete)
        for row in rows_to_delete:
            self.rows.remove(row)
        if self.eagerwrite:
            self.write()
        return len(rows_to_delete)

    def delete_row_func(self, column, func):
        if not self.initialized:
            self.init()

        rows_to_delete = []
        for row in self.rows:
            if func(row[column]):
                rows_to_delete.append(row)
#        map(self.rows.remove, rows_to_delete)
        for row in rows_to_delete:
            self.rows.remove(row)
        if self.eagerwrite:
            self.write()
        return len(rows_to_delete)

    def find_value(self, column, func):
        if not self.initialized:
            self.init()
        value = None
        for row in self.rows:
            if value is None:
                value = row[column]
                continue
            value = func(value, row[column])
        return value

    def find_row(self, column, func):
        if not self.initialized:
            self.init()
        value = None
        for row in self.rows:
            if value is None:
                value = row
                continue

            if func(row[column],value[column])==row[column]:
                value = row
        return value

    def min_value(self, column):
        return self.find_value(column, min)
    def min_row(self, column):
        return self.find_row(column, min)
    def max_value(self, column):
        return self.find_value(column, max)
    def max_row(self, column):
        return self.find_row(column, max)

    def get_row(self, column, value):
        if not self.initialized:
            self.init()

        for row in self.rows:
            if row[column] == value:
                return row

        return None

    def get_rows(self, column, value):
        if not self.initialized:
            self.init()

        result = []
        for row in self.rows:
            if row[column] == value:
                result.append(row)

        return result

    def get_column(self, column):
        if not self.initialized:
            self.init()
        results = []
        for row in self.rows:
            results.append(row[column])
        return results

if __name__=='__main__':
    import doctest
    doctest.testmod()
