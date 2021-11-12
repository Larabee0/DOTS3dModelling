using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Entities;

public struct EDConnectedTriangles : IBufferElementData
{
    public int ConnectedTriangleIndex;
    public Entity ConnectedTriangleEntity;
}