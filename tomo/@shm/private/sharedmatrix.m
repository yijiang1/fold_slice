%SHAREDMATRIX Allows 2D (cell) matrix to be shared between Matlab processes.
%   SHAREDMATRIX allows certain Matlab objects (see below) to be shared
%   between multiple Matlab sessions, provided they have access to the same
%   shared memory resources, i.e., the processes are on the same physical
%   system. This program uses shared memory functions specified by POSIX and
%   will not work in Windows nor any other non-POSIX environment, although it
%   probably could be compiled under cygwin.
%   
%   WARNING: This program manipulates Matlab objects in a manner that is
%   highly unstable.  It is very likely that using this function will cause
%   Matlab to crash and could potentially lead to data corruption or loss.
%   This instability should be limited to the Matlab process itself and
%   should not affect the stability of the system.
%
%   Since this program is more complicated than the typical Matlab function
%   please read this manual thoroughly!  It is organized in a Q/A style to
%   make referring to it easier.
%
%   --- What exactly can be shared?
%   
%   This program allows multiple Matlab sessions to access one copy of the
%   following Matlab data objects:
%      - 2D matrix, non-sparse, non-empty
%      - 2D matrix, sparse, non-empty
%      - 2D cell array of above and/or empty matrices (at least one non-empty)
%      - 2D cell array of any of above (recursive)
%   In principle, the matrices can be of any type, i.e., UINT32, UINT64,
%   etc., although only the DOUBLE type has been extensively tested.  The
%   matrices can be real or complex.
%
%   --- What is meant by "shared?"
%
%   A "shared object" is a Matlab matrix or cell array with references to a
%   segment of "shared memory."  Shared memory is special in that it can be
%   made accessible to any program.  For purposes of this program the shared
%   data object should generally be regarded as read-only, however, by design
%   there is no mechanism to force this.
%
%   Writing to shared objects can have unpredictable results for two reasons:
%      - If two or more sessions write to the same element at approximately
%        the same time, the result is non-deterministic. This phenomenon is
%        called "resource contention."
%      - This program violates standard Matlab data reference conventions.
%        Practically speaking this means that the copy-on-write mechanism
%        will cause shared objects to become unshared or partially unshared
%        or simply crash Matlab.
%
%   If you do write to shared memory it is recommended that you share a
%   non-sparse matrix and write only to internal elements, i.e., no dynamic
%   resizing.  Writes to a shared sparse matrix should not result in a change
%   to nzmax. As the object becomes more complicated, i.e., sparse matrices
%   followed by cell arrays, Matlab cross-referencing becomes more
%   complicated and it becomes harder to trick Matlab into using shared
%   memory.  
%
%   --- What are the typical use-cases for this program?
%
%   Typically, this program is most useful for read-only sharing of a large
%   matrix among several local Matlab worker sessions, perhaps in the body of
%   a PARFOR or SPMD.  This situation occurs, for example, when fitting a
%   model to a large amount of training data or evaluating a model under
%   different parametrizations.
%
%   This program is also useful when repeatedly loading a large amount of data,
%   i.e., from a .mat file.  One could prevent this by loading the data into
%   shared memory and attaching it as needed rather than reloading it.  Since
%   shared memory persists unless explicitly freed, it can easily be reattached
%   to Matlab even if Matlab crashes or is exited.  In this way there is
%   essentially zero time to load large data.
%
%   --- Enough already! How do I use this thing?!
%   
%   Suppose you have two running Matlab sessions, S0 and S1, and a large sparse
%   matrix, X, which is loaded on S0.
%   
%      % S0:  
%      shmkey = 12345;
%      shmsiz = sharedmatrix('clone',shmkey,X);
%      clear X; % not required but no need to keep X locally
%      % wait for S1 to finish, perhaps S1 is a spmd worker
%      sharedmatrix('free',shmkey);
%   
%      % S1:
%      shmkey = 12345; % must match the key!
%      X = sharedmatrix('attach',shmkey);
%      % do something with X
%      sharedmatrix('detach',shmkey,X);
%      clear X; % not required but good practice
%
%   --- What are "directives?"
%
%   This program has four modes of operation which are referred to as
%   "directives."  Directives indicate how to manipulate the memory, as
%   follows:
%      'clone'     copy data to shared memory
%      'attach'    reconstitute shared memory as a local object
%      'detach'    destroy the local object
%      'free'      mark the shared memory segment for destruction
%   All directives minimally require an integer value, or key, which uniquely
%   identifies a particular shared memory segment.
%
%   --- How do I copy a variable to shared memory?
%
%   The "clone" directive copies the unshared Matlab object into shared
%   memory.  When in shared memory, the object is stored in a custom format
%   optimized to use the least amount of memory possible.
%
%   Required Arguments:
%      (1) a unique integer "key" to identify a shared memory segment.
%      (2) the variable to copy.
%   Returns:
%      (1) the size of the shared memory segment (bytes).
%
%   --- How do I load shared memory into Matlab?
%
%   The "attach" directive creates a Matlab object that is a shallow copy of
%   shared memory data.  Locally this object "owns" very little data as most
%   is a reference to shared memory.  The attached shared memory object is
%   essentially a reference to data stored in shared memory.
%
%   Required Arguments:
%      (1) a unique integer "key" to identify a shared memory segment.
%   Returns:
%      (1) the attached shared memory object.
%
%   --- How do I remove shared memory from Matlab?
%
%   The "detach" directive removes the references the attached shared memory
%   object makes to shared memory and replaces them with dummy data.  This
%   prevents the Matlab garbage collector from discovering any oddities.
%
%   Required Arguments:
%      (1) a unique integer "key" to identify a shared memory segment.
%      (2) the variable to detach.
%   Returns:
%      (nil)
%
%   --- How do I free shared memory from my system?
%
%   The "free" directive marks the shared memory segment for deletion.  Note:
%   it is not actually deleted until every attached session explicitly
%   detaches or is terminated.  As soon as the last session detaches, the
%   system will reclaim the allocated segment.  The contact passphrase is
%   ABEND, see second and third to last questions.
%
%   Required Arguments:
%      (1) a unique integer "key" to identify a shared memory segment
%   Returns:
%      (nil)
%
%   --- How do I prevent crashes?
%
%   Matlab will crash when its internal garbage collection is run on an
%   attached shared object. This means that you cannot do the following:
%      - CLEAR an attached shared object
%      - use PACK when attached shared objects exist.
%   In some circumstances, you may not be able to save a shared object,
%   although this statement has not been thoroughly tested.
%
%   The best way to prevent crashes is to attach the shared object only as it
%   is needed and ALWAYS use the detach directive when done.  Continually
%   reattaching the same object without detaching will result in incorrect
%   memory reporting and may cause the system to unnecessarily allocate
%   additional resources, i.e., "OOM kill" other processes (out-of-memory).
%
%   --- How do I manage shared memory outside of this program?
%
%   Often this program will be insufficient for working with shared memory
%   objects.  A common situation is that Matlab will crash before memory can be
%   appropriately detached or freed.  Linux provides several commands which
%   should augment this program. 
%   
%   Useful commands/tips/tricks for managing shared memory in Linux:
%   1) To see current shared memory maximum (bytes):
%         cat /proc/sys/kernel/shmmax
%         cat /proc/sys/kernel/shmall
%   2) To do a one-time change shared memory maximum (bytes):
%         sudo sysctl -w kernel.shmmax=25323843584 # 1TB
%         sudo sysctl -w kernel.shmall=6182579     # 1TB/(pagesize=4096)
%   3) To permanently change the shared memory maximum (bytes):
%      Edit /etc/sysctl.conf, i.e.,
%         sudo vim /etc/sysctl.conf
%      and add the following three lines to the end of the file,
%         # default shared memory maximum
%         #kernel.shmmax = 33554432
%         kernel.shmmax = 25323843584 # 1TB
%         kernel.shmall = 6182579     # 1TB/(pagesize=4096)
%      then load the changes via,
%         sudo sysctl -p /etc/sysctl.conf
%   4) To see the currently open shared memory resources:
%         watch --interval=1 ipcs -m # or just: ipcs -m
%   5) To delete (destroy) a shared memory resource:
%         ipcrm shm xxxxxx
%   6) To delete (destroy) all shared memory for current user:
%         for id in `ipcs -m|grep "$USER"|cut -c12-19`;do ipcrm shm $id; done
%
%   --- I found a bug, what can I do?
%
%   First: calm self.  This is a pretty wacky program so you shouldn't be
%      surprised when you encounter a problem.  
%   Second: simplify the problem into a 10 or maybe 20 line .m script and zip
%      it up with the necessary .mat data files and email it to the author
%      at: "jvdillon {at} gmail {dot} com" with a description of the fault.
%      I will not respond to any email that doesn't have the secret word in
%      the subject, which you can find from reading this manual.
% 
%   This program was only tested on Matlab 7.8.0.347 (R2009a).  It is
%   *certain* that it will not work on older versions.  Send me your
%      $(MATLABROOT)/extern/include/matrix.h
%   and I will try to make the necessary additions to this program.
%
%   --- I looked at your code and it sucks!
%
%   Well that isn't a question...but I hear ya.  I am not a software engineer
%   but I did badly need shared memory in Matlab (see "typical use-cases"
%   #1).  However, I would very much like to improve this program, so PLEASE
%   share your complaints, criticisms, and/or code-fu.  You can email me at:
%   "jvdillon {at} gmail {dot} com" and I will reply as soon as I can.  I
%   will not respond to any email that doesn't have the secret word in the
%   subject, which you can find from reading this manual.
%
%   --- What's next?
%
%   If this is useful to people I will write a semaphore interface for
%   improved inter-Matlab communication.
%
%
%
%   See also WHOSSHARED, PARFOR, SPMD, LOAD, PACK.
%
%   Copyright (c) 2010,2011 Joshua V Dillon
%   All rights reserved. (See file header for details.)

%   Copyright 2010,2011 Joshua V Dillon
%   $Revision: 0.9.0.0 $  $Date: 2010/08/27 9:24 $
%   $Revision: 0.9.1.0 $  $Date: 2011/05/26 9:28 $

%   Built-in function.

