## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2024
## ...........................................................................

from .latticeExceptions import MoveSetException

class MoveSet:
    """
    MoveSet objects define a move made which can be used to improve perfomance when determining where
    a new structure is going to lie.

    Basically, this has no functionality but is a nice structured wrapper around some information 
    which is often useful. Just makes the code/logic much clearer than bungling everything into a 
    list or dictionary.
    """

    def __init__(self, moveType, dimensions, rotation_degrees=None, rotation_axis=None, rotation_anchor=None, translation_offset=None):
        """
        Initialize and validate a MoveSet description.

        MoveSet objects only contain variables - no functionality. All
        variables must be set on initialization, and the constructor validates
        that the combination of arguments is internally consistent for the
        requested ``moveType`` (raising ``MoveSetException`` otherwise).

        Parameters
        ----------
        moveType : str
            The type of move being described. Must be either ``'rotation'`` or
            ``'translation'``.

        dimensions : list or array-like
            Real dimensions of the system being moved (i.e. an n-dimensional
            array containing the lattice size along each axis).

        rotation_degrees : int, optional
            Number of degrees the rotation operation is carrying out. Must be
            one of 90, 180 or 270. Required for (and only used by) rotation
            moves.

        rotation_axis : {'x', 'y', 'z'}, optional
            Axis around which rotation is occurring. Only used (and required)
            for 3D rotation moves.

        rotation_anchor : list or array-like, optional
            Position on the grid used to offset positions back to the origin
            for rotation. This is usually set to the center of mass of the
            chain(s) being rotated. Defining it explicitly is useful because
            for cluster rotations a random assortment of positions from within
            the cluster and at its interface are used and the COM is only an
            approximation; the anchor ensures every position has exactly the
            same translation/rotation operation applied. Required for rotation
            moves and must match the dimensionality of ``dimensions``.

        translation_offset : list or array-like, optional
            n-dimensional array giving the translation offset being applied.
            Required for translation moves and must match the dimensionality of
            ``dimensions``.

        Returns
        -------
        None

        Raises
        ------
        MoveSetException
            If ``moveType`` is not ``'rotation'`` or ``'translation'``, or if
            the arguments required for the requested move are missing,
            invalid, or of the wrong dimensionality.
        """
        
        if moveType not in ['rotation', 'translation']:
            raise MoveSetException('Passed an invalid moveset')

        self.moveType           = moveType
        self.dimensions         = dimensions

        
        if moveType == 'translation':

            if translation_offset is None:
                raise MoveSetException('translation move requires a translation_offset')

            if not len(dimensions) == len(translation_offset):                
                raise MoveSetException('translation_offset must be the same number of dimensions as the dimensions variable')
                

        if moveType == 'rotation':

            # check we have a rotation_degree
            if rotation_degrees is None:
                raise MoveSetException('rotation move requires a rotation_degree variable (90,180,270)')

            # check it's valid
            if rotation_degrees not in [90, 180, 270]:
                raise MoveSetException('rotation_degree variable must be a cardinal rotation direction (90,180,270)')

            # check we have a rotation anchor point
            if rotation_anchor is None:
                raise MoveSetException('rotation_anchor required for a rotation move')

            # check rotation anchor position is the correct dimensionality
            if not len(rotation_anchor) == len(dimensions):
                raise MoveSetException('rotation_anchor must be the same dimensionality as the dimensions keyword')

            # if we're in 3D need a rotation axis
            if len(dimensions) == 3 and rotation_axis is None:
                raise MoveSetException('For 3D rotation must specificy a rotation_axis')
                
            # check rotation axis is a valid option
            if len(dimensions) == 3 and rotation_axis not in ['x', 'y', 'z']:
                raise MoveSetException("rotation_axis must be one of 'x', 'y' or 'z'")
                                                

        self.rotation_degrees   = rotation_degrees
        self.rotation_axis      = rotation_axis
        self.translation_offset = translation_offset
        self.rotation_anchor    = rotation_anchor
        
        

            
