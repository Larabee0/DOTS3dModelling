using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Entities;

public struct VDRelatedVert : IBufferElementData
{
    public int RelatedVert;
}

public struct VDConnectedVert : IBufferElementData
{
    public int ConnectedVertIndex;
    public Entity ConnectedVertEntity;
}

public struct VDConnectedEdge : IBufferElementData
{
    public int ConnectedEdgeIndex;
    public Entity ConnectedEdgeEntity;
}

public struct VDConnectedTriangle : IBufferElementData
{
    public int ConnectedTriangle;
    public Entity ConnectedTriangleEntity;
}
