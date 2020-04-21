Epitalnew
=======

Goal
-----------

This is the new version of my code to calculate quantum levels, wavefunctions and some properties of semiconductor heterostructures (mainly grown by epitaxy).

Initially it was focused on  III-V materials, expanded to II-VI. Some VdW materials support is planned.

The documentation started well but then I mislead. Sorry! I promise to improve it.

License
-----------

Please, give me some acknowledgement if you publish something which use this code.


Requirements
-----------

I recommend to use CMake to compile it (if you are using VTK it is almost mandatory).

It requires armadillo and Eigen3 for some linear algebra and matrix operations.
It will be nice to have OpenMP since it has multi-thread support through it, and most of the operations scale with the number of processors.
VTK is required for some 3D operations (Heterostructure3D.hpp). Sometimes is trick to install it. If you have problems quit Heterostructure3D.hpp/Heterostructure3D.hpp

How to use
-----------

The folder /src/Examples has some examples of tested simulations


To compile:

1. Create a folder inside /src with a name .build
2. Run CMake and check if all the dependencies are satisfied
3. Run make and Good Luck!

Example:

```
cd src
mk .build
cd .build
cmake ..
build
```

TODO
-----------
1. Update Docs
2. Do test routines
3. Improve examples
4. VdW materials InSe, GaSe and MoSe2






Ignore it:
-----------
.md In-situ Reference ;)


Two spaces at the end of a line  
produces a line break.

Text attributes _italic_,
**bold**, `monospace`.

Horizontal rule:

---

Strikethrough:
~~strikethrough~~

Bullet list:

  * apples
  * oranges
  * pears

Numbered list:

  1. lather
  2. rinse
  3. repeat

An link[Google](http://google.com).

![Image](Icon-pictures.png "icon")

> Markdown uses email-style > characters for blockquoting.

Inline <abbr title="Hypertext Markup Language">HTML</abbr> is supported.
