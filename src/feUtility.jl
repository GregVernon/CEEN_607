using NamedDims
import LinearAlgebra

function affineTransformationMatrix_2D(;translate="TranslateVector",scale="ScaleVector",rotate="RotateAngle",reflect="Reflection",shear="ShearVector")
    ReferenceMatrix = NamedDimsArray{(:ℝᴺ,:ℝᴺ)}(float([1 0 0; 0 1 0; 0 0 1]))
    TranslateMatrix = NamedDimsArray{(:ℝᴺ,:ℝᴺ)}(float([1 0 0; 0 1 0; 0 0 1]))
    ScaleMatrix = NamedDimsArray{(:ℝᴺ,:ℝᴺ)}(float([1 0 0; 0 1 0; 0 0 1]))
    RotateMatrix = NamedDimsArray{(:ℝᴺ,:ℝᴺ)}(float([1 0 0; 0 1 0; 0 0 1]))
    ReflectMatrix = NamedDimsArray{(:ℝᴺ,:ℝᴺ)}(float([1 0 0; 0 1 0; 0 0 1]))
    ShearMatrix = NamedDimsArray{(:ℝᴺ,:ℝᴺ)}(float([01 0 0; 0 1 0; 0 0 1]))

    # Translate Matrix
    if typeof(translate) == String
        nothing
    else
        TranslateMatrix[1,3] = translate[1]
        TranslateMatrix[2,3] = translate[2]
    end

    # Scale Matrix
    if typeof(scale) == String
        nothing
    else
        ScaleMatrix[1,1] = scale[1] - 1.0
        ScaleMatrix[2,2] = scale[2] - 1.0
    end

    # Rotation Matrix
    if typeof(rotate) == String
        nothing
    else
        rotate = -1.0 * rotate
        RotateMatrix[[1,2],[1,2]] .= [cos(rotate) sin(rotate); -sin(rotate) cos(rotate)]
    end

    # Reflection Matrix
    if typeof(reflect) == String
        nothing
    else
        ReflectMatrix[1,1] = 1.0 - 2.0*reflect[1]
        ReflectMatrix[2,2] = 1.0 - 2.0*reflect[2]
    end

    # Shear Matrix
    if typeof(shear) == String
        nothing
    else
        ShearMatrix[1,2] = tan(shear[1])
        ShearMatrix[2,1] = tan(shear[2])
    end

    # Combine
    TransformMatrix = (ReferenceMatrix * TranslateMatrix * ScaleMatrix * RotateMatrix * ReflectMatrix * ShearMatrix) 
    return TransformMatrix
end