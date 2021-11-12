using UnityEngine;
using Unity.Mathematics;
using Unity.Entities;
using Unity.Collections;
using System;

public struct MeshManagerType : IComponentData { }
public struct MeshManagerPreStart : IComponentData { }
public struct MeshManagerStart : IComponentData { }
public struct MeshManagerStartVert : IComponentData { }
public struct MeshManagerStartEdge : IComponentData { }
public struct MeshManagerStartTri : IComponentData { }
public struct MeshManagerStartedVert : IComponentData { }
public struct MeshManagerStartedEdge : IComponentData { }
public struct MeshManagerStartedTri : IComponentData { }
public struct MeshManagerStarted : IComponentData { }
public struct UpdateVisualComponent : IComponentData { }
public struct Vertices : IBufferElementData
{
    public Vector3 GlobalVertex;
    public Vector3 LocalVertex;
    public int SimplifiedVertex;
}

public struct Triangles : IBufferElementData
{
    public int GlobalTriangle;
    public int SimplifiedTriangle;
}
public struct SimplifiedVertices : IBufferElementData
{
    public Entity VertEntity;
    public Vector3 Vertex;
}

public struct SimplifiedEdges : IBufferElementData
{
    public int EdgeIndex;
    public int SimpVert1Index;
    public int SimpVert2Index;
    public Entity SimpVert1Entity;
    public Entity SimpVert2Entity;
    public Entity EdgeEntity;
    public float3 EdgePosition;
}
public struct SimplifiedTriangles : IBufferElementData
{
    public int TriangleIndex;
    public int SimpVert1Index;
    public int SimpVert2Index;
    public int SimpVert3Index;
    public Entity SimpVert1Entity;
    public Entity SimpVert2Entity;
    public Entity SimpVert3Entity;
    public Entity TriangleEntity;
}