using Unity.Entities;
using Unity.Mathematics;

public struct TriangleDataStore : IComponentData
{
    public int TriangleIndex;
    public float3 TrianglePosition;
}
