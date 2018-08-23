# Memory Management in pawpyseed (Don't worry, it's not scary!)

The goal of pawpyseed is to combine power and flexibility for numerical
Kohn-Sham state analysis. To this end, pawpyseed implements a C backend
for heavy parallel calculations. In addition, due to the large amount of
data in plane-wave basis set wavefunctions, much of the data used
is stored in the C backend to avoid unnecessary C to Python data passing.

One downside of this system is that pawpyseed requires a little bit
of user memory management. If you only handle a couple WAVECAR files
per script, you probably don't even need to worry about it, but
for scripts that process many WAVECAR files
or very large WAVECAR files, you might run out of memory
if you aren't careful! This tutorial shows how to manage memory in pawpyseed.
It's short because there isn't much to know: you just have to call
Wavefunction.free_all() and/or Projector.free_all() at the right time.

## The Wavefunction.free_all() method

Say you made a wavefunction:

```
wf = Wavefunction.from_directory('.')
```

As is, pawpyseed's C backend contains two pointers. One is for the variable `projector_list`,
which contains all the data about the projectors and partial waves for the
wavefunction. The other is `pwf_ptr`, which contains the pseudowavefunction
plane-wave constants. When you are done using `wf` (**and not before then**),
call

```
wf.free_all()
```

to take care of the data associated with those two pointers. After this, `wf` is no
longer usable, and you don't have any memory leaks!

## The Projector.free_all() method

The Projector.free_all() method is a little more nuanced, but still easy to use.
Let's say you have

```
wf = Wavefunction.from_directory('defect', setup_projectors=False)
basis = Wavefunction.from_directory('bulk', setup_projectors=False)
pr = Projector(wf, basis)
```

(See Tutorial 2 for why you want to set setup_projectors=False).
Let's say you then make use of the `Projector` object you set up, for example:

```
res = pr.defect_band_analysis()
```

And now you're done with it. **Projector.free_all() only frees the combo projector_list
used by the Projector object. It does not free pseudowavefunctions or projector lists
used by the constituent Wavefunction objects**. However, the wavefunction objects might
rely on the `projector_list` stored in the `Projector` object, depending on how
they were set up. Therefore, **always call free_all() for Wavefunction objects before
the associated Projector object!**.

**GOOD**
```
basis.free_all()
wf.free_all()

pr.free_all()
```

**BAD**
```
pr.free_all()

basis.free_all()
wf.free_all()
```

