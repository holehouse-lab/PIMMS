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
        
        MoveSet objects only contain variables - no functionality. Variables must be set on initialization


        moveType:
        type of move being described ('rotation' or 'translation')
        
        dimensions:
        Real dimensions of system being moved (i.e. a n-dimensional array containing lattice size)

        rotation_degrees
        Number of degrees a rotation operation is carrying out (90/180/270)

        rotation_axis:
        axis around which rotation is occuring (only used for 3D rotation) - (x/y/z)

        rotation_anchor:
        position on grid used to offset the positions to re-set them to the origin for rotation - this is
        usually set as the center of mass of the chain(s) being rotated, but defining it here is useful because
        the COM is an approximation so when we carry out cluster rotations later and use a random assortment of positions 
        from within the cluster and at the clutser interface the COM will change. A rotation_anchor ensures all
        rotation operations apply exactly the same set of translation/rotation operations.

        translation_offset:
        n-dim array giving the translation offset being used



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
        
        

            
