using Unity.Entities;

public struct TDConnectedEdges : IBufferElementData
{
    public int ConnectedEdgeIndex;
    public Entity ConnectedEdgeEntity;
}
