using System.Collections;
using System.Collections.Generic;
using Unity.Entities;
using UnityEngine;
using Unity.Transforms;
using Unity.Mathematics;
using Unity.Collections;
using Unity.Collections.LowLevel.Unsafe;
using Unity.Rendering;
using Unity.Jobs;
using Unity.Burst;
using Unity.Physics;
using System.Linq;
#region UnusedSystems
[DisableAutoCreation]
/// <summary>
/// this is the best we can do for none burst
/// </summary>
public class VisualsUpdate : ComponentSystem
{
    protected override void OnUpdate()
    {
        EntityQuery MeshManager = GetEntityQuery(typeof(MeshManagerType));
        NativeArray<Entity> MeshManagerEntity = MeshManager.ToEntityArray(Allocator.Temp);

        BufferFromEntity<SimplifiedVertices> SimpVertFromEntity = GetBufferFromEntity<SimplifiedVertices>();
        BufferFromEntity<SimplifiedEdges> SimpEdgeFromEntity = GetBufferFromEntity<SimplifiedEdges>();

        DynamicBuffer<SimplifiedVertices> simplifiedVertices = SimpVertFromEntity[MeshManagerEntity[0]];
        DynamicBuffer<SimplifiedEdges> simplifiedEdges = SimpEdgeFromEntity[MeshManagerEntity[0]];
        MeshManagerEntity.Dispose();
        Entities.WithAll<EdgeDataStore>().ForEach((ref EdgeDataStore edgeDataStore, ref Translation translation, ref Rotation rotation, ref NonUniformScale scale) =>
        {
            SimplifiedEdges simplifiedEdge = simplifiedEdges[edgeDataStore.EdgeIndex];
            float3 A = simplifiedVertices[simplifiedEdge.SimpVert1Index].Vertex;
            float3 B = simplifiedVertices[simplifiedEdge.SimpVert2Index].Vertex;
            float3 EdgeAveragePositon = ((A + B) / 2);
            //float3 C = A - B;
            translation.Value = EdgeAveragePositon;
            scale.Value = new float3(0.025f, 0.025f, math.distance(A,B));
            rotation.Value = quaternion.LookRotationSafe(A - EdgeAveragePositon, new float3(1, 1, 1));
        });


        Entities.WithAll<VertDataStore>().ForEach((ref Translation translation, ref VertDataStore VDS) =>
        {
            translation.Value = new float3(simplifiedVertices[VDS.SimplifiedVertexIndex].Vertex);
        });
    }
}
[DisableAutoCreation]
/// <summary>
/// This is terrible
/// </summary>
public class UpdateVisualsSystemJobChunk : SystemBase
{
    private EndSimulationEntityCommandBufferSystem endSimulationEntityCommandBufferSystem;
    private EntityQueryDesc MeshManagerQuery;
    private EntityQueryDesc VertDataStoreQuery;
    private EntityQueryDesc EdgeDataStoreQuery;
    protected override void OnCreate()
    {
        endSimulationEntityCommandBufferSystem = World.GetOrCreateSystem<EndSimulationEntityCommandBufferSystem>();
        MeshManagerQuery = new EntityQueryDesc { All = new ComponentType[] { typeof(MeshManagerType) } };
        VertDataStoreQuery = new EntityQueryDesc { All = new ComponentType[] { typeof(VertDataStore), typeof(Translation) } };
        EdgeDataStoreQuery = new EntityQueryDesc { All = new ComponentType[] { typeof(EdgeDataStore), typeof(Translation), typeof(Rotation), typeof(NonUniformScale) } };
    }

    protected override void OnUpdate()
    {
        BufferFromEntity<SimplifiedVertices> SimpVertFromEntity = GetBufferFromEntity<SimplifiedVertices>();
        BufferFromEntity<SimplifiedEdges> SimpEdgeFromEntity = GetBufferFromEntity<SimplifiedEdges>();
        EntityQuery MeshManager = GetEntityQuery(MeshManagerQuery);
        NativeArray<Entity> MeshManagerEntity = MeshManager.ToEntityArray(Allocator.TempJob);
        DynamicBuffer<SimplifiedVertices> simplifiedVertices = SimpVertFromEntity[MeshManagerEntity[0]];
        DynamicBuffer<SimplifiedEdges> simplifiedEdges = SimpEdgeFromEntity[MeshManagerEntity[0]];
        MeshManagerEntity.Dispose();

        EntityQuery VertDataStore = GetEntityQuery(VertDataStoreQuery);
        NativeArray<Entity> VertEntities = VertDataStore.ToEntityArray(Allocator.TempJob);
        NativeArray<VertDataStore> targetVDSArray = VertDataStore.ToComponentDataArray<VertDataStore>(Allocator.TempJob);
        for (int i = 0; i < VertEntities.Length; i++)
        {
            endSimulationEntityCommandBufferSystem.CreateCommandBuffer().SetComponent<Translation>(VertEntities[i], new Translation { Value = simplifiedVertices[targetVDSArray[i].SimplifiedVertexIndex].Vertex });
        }
        VertEntities.Dispose();
        targetVDSArray.Dispose();


        EntityQuery EdgeDataStore = GetEntityQuery(EdgeDataStoreQuery);
        int arraylengths = EdgeDataStore.CalculateEntityCount();
        NativeArray<Entity> EdgeEntities = EdgeDataStore.ToEntityArray(Allocator.TempJob);
        
        NativeArray<EdgeDataStore> targetEDSArray = EdgeDataStore.ToComponentDataArray<EdgeDataStore>(Allocator.TempJob);
        NativeArray<float3> EDSposArray = new NativeArray<float3>(arraylengths, Allocator.TempJob);
        NativeArray<float3> EDSscaleArray = new NativeArray<float3>(arraylengths, Allocator.TempJob);
        NativeArray<quaternion> EDSrotArray = new NativeArray<quaternion>(arraylengths, Allocator.TempJob);

        UpdateEdgeVisual job = new UpdateEdgeVisual
        {
            EDS = targetEDSArray,
            SimplifiedVertices = simplifiedVertices,
            SimplifiedEdges = simplifiedEdges,
            Position = EDSposArray,
            Scale = EDSscaleArray,
            Rotation = EDSrotArray
        };
        JobHandle handle = job.Schedule(EdgeDataStore);
        handle.Complete();
        targetEDSArray.Dispose();
        for (int i = 0; i < EdgeEntities.Length; i++)
        {
            endSimulationEntityCommandBufferSystem.CreateCommandBuffer().SetComponent<Translation>(EdgeEntities[i], new Translation { Value = EDSposArray[i] });
            endSimulationEntityCommandBufferSystem.CreateCommandBuffer().SetComponent<NonUniformScale>(EdgeEntities[i], new NonUniformScale { Value = EDSscaleArray[i] });
            endSimulationEntityCommandBufferSystem.CreateCommandBuffer().SetComponent<Rotation>(EdgeEntities[i], new Rotation { Value = EDSrotArray[i] });
        }
        EdgeEntities.Dispose();
        EDSposArray.Dispose();
        EDSscaleArray.Dispose();
        EDSrotArray.Dispose();
    }
    [BurstCompile]
    private struct UpdateEdgeVisual : IJobChunk
    {
        [ReadOnly]
        public NativeArray<EdgeDataStore> EDS;
        [ReadOnly]
        public DynamicBuffer<SimplifiedVertices> SimplifiedVertices;
        [ReadOnly]
        public DynamicBuffer<SimplifiedEdges> SimplifiedEdges;

        public NativeArray<float3> Position;
        public NativeArray<float3> Scale;
        public NativeArray<quaternion> Rotation;
        public void Execute(ArchetypeChunk chunk, int chunkIndex, int firstEntityIndex)
        {
            for (int i = 0; i < chunk.Count; i++)
            {
                SimplifiedEdges simplifiedEdge = SimplifiedEdges[EDS[i].EdgeIndex];
                float3 A = SimplifiedVertices[simplifiedEdge.SimpVert1Index].Vertex;
                float3 B = SimplifiedVertices[simplifiedEdge.SimpVert2Index].Vertex;
                float3 EdgeAveragePositon = (A + B) / 2;
                Position[i] = EdgeAveragePositon;
                Scale[i] = new float3(0.025f, 0.025f, math.distance(A, B));
                Rotation[i] = quaternion.LookRotationSafe(A - EdgeAveragePositon, new float3(1, 1, 1));
            }
        }
    }
}
[DisableAutoCreation]
/// <summary>
/// This works really well with burst enabled and the same as VisualUpdate when burst is disabled.
/// </summary>
public class UpdateVisualsSystemJobEntityBatch : SystemBase
{
    private EntityQueryDesc MeshManagerQuery;
    private EntityQueryDesc VertDataStoreQuery;
    private EntityQueryDesc EdgeDataStoreQuery;
    protected override void OnCreate()
    {
        MeshManagerQuery = new EntityQueryDesc { All = new ComponentType[] { typeof(MeshManagerType) } };
        VertDataStoreQuery = new EntityQueryDesc { All = new ComponentType[] { typeof(VertDataStore), typeof(Translation) } };
        EdgeDataStoreQuery = new EntityQueryDesc { All = new ComponentType[] { typeof(EdgeDataStore), typeof(Translation), typeof(Rotation), typeof(NonUniformScale) } };
    }

    protected override void OnUpdate()
    {
        BufferFromEntity<SimplifiedVertices> SimpVertFromEntity = GetBufferFromEntity<SimplifiedVertices>();
        BufferFromEntity<SimplifiedEdges> SimpEdgeFromEntity = GetBufferFromEntity<SimplifiedEdges>();
        EntityQuery MeshManager = GetEntityQuery(MeshManagerQuery);
        NativeArray<Entity> MeshManagerEntity = MeshManager.ToEntityArray(Allocator.TempJob);
        DynamicBuffer<SimplifiedVertices> simplifiedVertices = SimpVertFromEntity[MeshManagerEntity[0]];
        DynamicBuffer<SimplifiedEdges> simplifiedEdges = SimpEdgeFromEntity[MeshManagerEntity[0]];
        MeshManagerEntity.Dispose();

        UpdateVertexVisuals vertexJob = new UpdateVertexVisuals
        {
            VDSTypeHandle = this.GetComponentTypeHandle<VertDataStore>(true),
            translationTypeHandle = this.GetComponentTypeHandle<Translation>(),
            SimplifiedVertices = simplifiedVertices
        };
        EntityQuery VertQuery = GetEntityQuery(VertDataStoreQuery);
        JobHandle vertexHandle = vertexJob.Schedule(VertQuery, this.Dependency);
        //vertexHandle.Complete();
        EntityQuery EdgeQuery = GetEntityQuery(EdgeDataStoreQuery);
        UpdateEdgeVisuals edgeJob = new UpdateEdgeVisuals
        {
            EDSTypeHandle = this.GetComponentTypeHandle<EdgeDataStore>(true),
            translationTypeHandle = this.GetComponentTypeHandle<Translation>(),
            scaleTypeHandle = this.GetComponentTypeHandle<NonUniformScale>(),
            RotationTypeHandle = this.GetComponentTypeHandle<Rotation>(),
            SimplifiedVertices = simplifiedVertices,
            SimplifiedEdges = simplifiedEdges
        };
        JobHandle edgeHandle = edgeJob.Schedule(EdgeQuery, vertexHandle);
        
        edgeHandle.Complete();
    }


}
#endregion

#region StartUpSystem

public class PreAndPostStartUpSystem : ComponentSystem
{
    private EntityManager entityManager;
    private EntityArchetype VertDataStoreArch;
    private EntityArchetype TriangleDataStoreArch;
    private EntityArchetype EdgeDataStoreArch;
    private EntityCommandBuffer StartUpCommandBuffer;
    private EndSimulationEntityCommandBufferSystem endSimulationEntityCommandBufferSystem;
    
    private EntityQueryDesc TriangleDataStoreQuery;
    private EntityQueryDesc EdgeDataStoreQuery;
    private EntityQueryDesc VertDataStoreQuery;
    protected override void OnCreate()
    {
        endSimulationEntityCommandBufferSystem = World.GetOrCreateSystem<EndSimulationEntityCommandBufferSystem>();
        
        //this.Enabled = false;
        entityManager = World.DefaultGameObjectInjectionWorld.EntityManager;
        VertDataStoreQuery = new EntityQueryDesc
        {
            All = new ComponentType[] { typeof(VertDataStore) }
        };
        EdgeDataStoreQuery = new EntityQueryDesc
        {
            All = new ComponentType[] { typeof(EdgeDataStore) }
        };
        TriangleDataStoreQuery = new EntityQueryDesc
        {
            All = new ComponentType[] { typeof(TriangleDataStore) }
        };
        
        VertDataStoreArch = entityManager.CreateArchetype(
            typeof(Translation),
            typeof(Scale), //
            typeof(RenderMesh),
            typeof(LocalToWorld),
            typeof(LocalToParent),//
            typeof(Parent),//
            typeof(RenderBounds),
            //typeof(PhysicsCollider),
            typeof(VertDataStore)
        ); 
        EdgeDataStoreArch = entityManager.CreateArchetype(
            typeof(Translation),
            typeof(NonUniformScale),
            typeof(Rotation),
            typeof(RenderMesh),
            typeof(LocalToWorld),
            typeof(LocalToParent),
            typeof(Parent),
            typeof(RenderBounds),
            //typeof(PhysicsCollider),
            typeof(EdgeDataStore)
        );

        TriangleDataStoreArch = entityManager.CreateArchetype(
            typeof(Translation),
            typeof(Rotation),
            typeof(LocalToWorld),
            typeof(LocalToParent),
            typeof(Parent),
            //typeof(PhysicsCollider),
            typeof(TriangleDataStore)
        );
    }
    public struct SimpVertContainer
    {
        public NativeList<int> RelatedVerts { get; set; }
        public int Index { get; set; }
        public SimpVertContainer(NativeList<int> relatedVerts, int index)
        {
            RelatedVerts = relatedVerts;
            Index = index;
        }
        
    }
    protected override void OnUpdate()
    {
        StartUpCommandBuffer = endSimulationEntityCommandBufferSystem.CreateCommandBuffer();

        Entities.WithAll<MeshManagerType, MeshManagerStartedVert, MeshManagerStartedEdge, MeshManagerStartedTri>().ForEach((Entity MeshManagerReal) =>
        {
            float SetUpFirst = UnityEngine.Time.realtimeSinceStartup;
            DynamicBuffer<SimplifiedVertices> SimplifiedVerticesBuffer = entityManager.GetBuffer<SimplifiedVertices>(MeshManagerReal);
            DynamicBuffer<SimplifiedEdges> SimplifiedEdgesBuffer = entityManager.GetBuffer<SimplifiedEdges>(MeshManagerReal);
            DynamicBuffer<SimplifiedTriangles> SimplifiedTrianglesBuffer = entityManager.GetBuffer<SimplifiedTriangles>(MeshManagerReal);
            NativeArray<SimplifiedVertices> SimplifiedVerticesArray = SimplifiedVerticesBuffer.ToNativeArray(Allocator.Temp);
            NativeArray<SimplifiedEdges> SimplifiedEdgesArray = SimplifiedEdgesBuffer.ToNativeArray(Allocator.Temp);
            NativeArray<SimplifiedTriangles> SimplifiedTrianglesArray = SimplifiedTrianglesBuffer.ToNativeArray(Allocator.Temp);
            //ConnectedEdges for vertices
            for (int i = 0; i < SimplifiedEdgesArray.Length; i++)
            {
                Entity Vert1Entitiy = SimplifiedVerticesBuffer[SimplifiedEdgesArray[i].SimpVert1Index].VertEntity;
                Entity Vert2Entitiy = SimplifiedVerticesBuffer[SimplifiedEdgesArray[i].SimpVert2Index].VertEntity;
                entityManager.GetBuffer<VDConnectedEdge>(Vert1Entitiy).Add(new VDConnectedEdge { ConnectedEdgeIndex = i, ConnectedEdgeEntity = SimplifiedEdgesArray[i].EdgeEntity });
                entityManager.GetBuffer<VDConnectedEdge>(Vert2Entitiy).Add(new VDConnectedEdge { ConnectedEdgeIndex = i, ConnectedEdgeEntity = SimplifiedEdgesArray[i].EdgeEntity });
                DynamicBuffer<VDConnectedVert> TempBuffer1 = entityManager.GetBuffer<VDConnectedVert>(Vert1Entitiy);
                DynamicBuffer<VDConnectedVert> TempBuffer2 = entityManager.GetBuffer<VDConnectedVert>(Vert2Entitiy);
                NativeArray<VDConnectedVert> VDConnectedVerts1Array = TempBuffer1.ToNativeArray(Allocator.Temp);
                NativeArray<VDConnectedVert> VDConnectedVerts2Array = TempBuffer2.ToNativeArray(Allocator.Temp);
                if (!VDConnectedVerts1Array.Contains(new VDConnectedVert { ConnectedVertIndex = SimplifiedEdgesArray[i].SimpVert2Index, ConnectedVertEntity = SimplifiedEdgesArray[i].SimpVert2Entity }))
                {
                    TempBuffer1.Add(new VDConnectedVert { ConnectedVertIndex = SimplifiedEdgesArray[i].SimpVert2Index, ConnectedVertEntity = SimplifiedEdgesArray[i].SimpVert2Entity });
                }
                if (!VDConnectedVerts2Array.Contains(new VDConnectedVert { ConnectedVertIndex = SimplifiedEdgesArray[i].SimpVert1Index, ConnectedVertEntity = SimplifiedEdgesArray[i].SimpVert1Entity }))
                {
                    TempBuffer2.Add(new VDConnectedVert { ConnectedVertIndex = SimplifiedEdgesArray[i].SimpVert1Index, ConnectedVertEntity = SimplifiedEdgesArray[i].SimpVert1Entity });
                }
                VDConnectedVerts1Array.Dispose();
                VDConnectedVerts2Array.Dispose();
            }
            //Connected vertices is already known.

            //ConnectedEdgesForTriangles
            for (int i = 0; i < SimplifiedTrianglesArray.Length; i++)
            {
                Entity Vert1Entitiy = SimplifiedVerticesBuffer[SimplifiedTrianglesArray[i].SimpVert1Index].VertEntity;
                Entity Vert2Entitiy = SimplifiedVerticesBuffer[SimplifiedTrianglesArray[i].SimpVert2Index].VertEntity;
                Entity Vert3Entitiy = SimplifiedVerticesBuffer[SimplifiedTrianglesArray[i].SimpVert3Index].VertEntity;

                entityManager.GetBuffer<VDConnectedTriangle>(Vert1Entitiy).Add(new VDConnectedTriangle { ConnectedTriangle = i, ConnectedTriangleEntity = SimplifiedTrianglesArray[i].TriangleEntity });
                entityManager.GetBuffer<VDConnectedTriangle>(Vert2Entitiy).Add(new VDConnectedTriangle { ConnectedTriangle = i, ConnectedTriangleEntity = SimplifiedTrianglesArray[i].TriangleEntity });
                entityManager.GetBuffer<VDConnectedTriangle>(Vert3Entitiy).Add(new VDConnectedTriangle { ConnectedTriangle = i, ConnectedTriangleEntity = SimplifiedTrianglesArray[i].TriangleEntity });

                DynamicBuffer<VDConnectedEdge> ConnectedEdgesVertexBuffer1 = entityManager.GetBuffer<VDConnectedEdge>(Vert1Entitiy);
                DynamicBuffer<VDConnectedEdge> ConnectedEdgesVertexBuffer2 = entityManager.GetBuffer<VDConnectedEdge>(Vert2Entitiy);
                DynamicBuffer<VDConnectedEdge> ConnectedEdgesVertexBuffer3 = entityManager.GetBuffer<VDConnectedEdge>(Vert3Entitiy);

                NativeArray<VDConnectedEdge> VDConnectedEdges1Array = ConnectedEdgesVertexBuffer1.ToNativeArray(Allocator.Temp);
                NativeArray<VDConnectedEdge> VDConnectedEdges2Array = ConnectedEdgesVertexBuffer2.ToNativeArray(Allocator.Temp);
                NativeArray<VDConnectedEdge> VDConnectedEdges3Array = ConnectedEdgesVertexBuffer3.ToNativeArray(Allocator.Temp);
                VDConnectedEdge[][] ConnectedEdges = new VDConnectedEdge[][] { VDConnectedEdges1Array.ToArray(), VDConnectedEdges2Array.ToArray(), VDConnectedEdges3Array.ToArray() };
                DynamicBuffer<TDConnectedEdges> TDConnectedEdges = entityManager.GetBuffer<TDConnectedEdges>(SimplifiedTrianglesArray[i].TriangleEntity);

                int secondaryK = 1;
                for (int k = 0; k < 3; k++)
                {
                    if (secondaryK == 3)
                    {
                        secondaryK = 0;
                    }
                    for (int h = 0; h < ConnectedEdges[k].Length; h++)
                    {
                        for (int l = 0; l < ConnectedEdges[secondaryK].Length; l++)
                        {
                            if (ConnectedEdges[k][h].ConnectedEdgeIndex == ConnectedEdges[secondaryK][l].ConnectedEdgeIndex)
                            {
                                TDConnectedEdges.Add(new TDConnectedEdges { ConnectedEdgeIndex = ConnectedEdges[k][h].ConnectedEdgeIndex, ConnectedEdgeEntity = SimplifiedEdgesArray[ConnectedEdges[k][h].ConnectedEdgeIndex].EdgeEntity });
                                DynamicBuffer<EDConnectedTriangles> EDConnectedTrinalgesBuffer = entityManager.GetBuffer<EDConnectedTriangles>(SimplifiedEdgesArray[ConnectedEdges[k][h].ConnectedEdgeIndex].EdgeEntity);
                                EDConnectedTrinalgesBuffer.Add(new EDConnectedTriangles { ConnectedTriangleIndex = i, ConnectedTriangleEntity = SimplifiedTrianglesArray[i].TriangleEntity });
                                goto Next;
                            }
                        }
                    }
                    Next:
                    secondaryK++;
                }
                VDConnectedEdges1Array.Dispose();
                VDConnectedEdges2Array.Dispose();
                VDConnectedEdges3Array.Dispose();
            }


            SimplifiedVerticesArray.Dispose();
            SimplifiedEdgesArray.Dispose();
            SimplifiedTrianglesArray.Dispose();
            StartUpCommandBuffer.RemoveComponent(MeshManagerReal, new ComponentTypes(typeof(MeshManagerStartedVert), typeof(MeshManagerStartedEdge), typeof(MeshManagerStartedTri)));
            StartUpCommandBuffer.AddComponent(MeshManagerReal, typeof(MeshManagerStarted));
            Debug.Log("PostStartData() Execution Time: " + (UnityEngine.Time.realtimeSinceStartup - SetUpFirst) * 1000f + "ms");
        });
        Entities.WithAll<MeshManagerType, MeshManagerPreStart>().ForEach((Entity MeshManagerReal) =>
        {
            float SetUpFirst = UnityEngine.Time.realtimeSinceStartup;
            #region Vertex Data Preperation
            DynamicBuffer<Vertices> VerticeBuffer = entityManager.GetBuffer<Vertices>(MeshManagerReal);
            DynamicBuffer<SimplifiedVertices> SimplifiedVerticesBuffer = entityManager.GetBuffer<SimplifiedVertices>(MeshManagerReal);
            NativeHashSet<Vector3> LocalVertsHashSet = new NativeHashSet<Vector3>(SimplifiedVerticesBuffer.Length, Allocator.Temp);
            foreach (Vertices vertex in VerticeBuffer)
            {
                LocalVertsHashSet.Add(vertex.LocalVertex);
            }
            NativeArray<Vector3> SimpVertTemp = LocalVertsHashSet.ToNativeArray(Allocator.Temp);
            LocalVertsHashSet.Dispose();

            int SimpVertTempLength = SimplifiedVerticesBuffer.Capacity = SimpVertTemp.Length;
            Dictionary<Vector3, List<VDRelatedVert>> relatedVertsDict = new Dictionary<Vector3, List<VDRelatedVert>>(SimpVertTempLength);
            for (int i = 0; i < SimpVertTemp.Length; i++)
            {
                SimplifiedVerticesBuffer.Add(new SimplifiedVertices { Vertex = SimpVertTemp[i] });
                relatedVertsDict.Add(SimpVertTemp[i], new List<VDRelatedVert>());
            }
            SimpVertTemp.Dispose();

            Dictionary<Vector3, SimpVertContainer> RelatedVertDict = new Dictionary<Vector3, SimpVertContainer>(SimpVertTempLength);
            for (int i = 0; i < SimplifiedVerticesBuffer.Length; i++)
            {
                RelatedVertDict.Add(SimplifiedVerticesBuffer[i].Vertex, new SimpVertContainer(new NativeList<int>(Allocator.Temp), i));
            }

            int VertBufferLength = VerticeBuffer.Length;
            NativeArray<int> VerticesToSimplifiedVerts = new NativeArray<int>(VertBufferLength, Allocator.Temp);
            Vector3 Temp;
            for (int i = 0; i < VertBufferLength; i++)
            {
                Temp = VerticeBuffer[i].LocalVertex;
                RelatedVertDict[Temp].RelatedVerts.Add(i);
                VerticesToSimplifiedVerts[i] = RelatedVertDict[Temp].Index;
                relatedVertsDict[Temp].Add(new VDRelatedVert { RelatedVert = i });
            }

            SimpVertContainer[] RelatedVertDictArray = RelatedVertDict.Values.ToArray();
            for (int i = 0; i < RelatedVertDictArray.Length; i++)
            {
                for (int k = 0; k < RelatedVertDictArray[i].RelatedVerts.Length; k++)
                {

                    Vertices aVertex = VerticeBuffer[k];
                    aVertex.SimplifiedVertex = RelatedVertDictArray[i].Index;
                    VerticeBuffer[k] = aVertex;
                }
            }
            #endregion

            #region Triangle Data Preperation
            DynamicBuffer<Triangles> TrianglesBuffer = entityManager.GetBuffer<Triangles>(MeshManagerReal);
            int length = TrianglesBuffer.Length;
            DynamicBuffer<SimplifiedTriangles> SimplifiedTrianglesBuffer = entityManager.GetBuffer<SimplifiedTriangles>(MeshManagerReal);
            int Triangles = 0;
            for (int i = 0; i < length; i += 3)
            {

                Triangles aTriangle = TrianglesBuffer[i];
                aTriangle.SimplifiedTriangle = VerticesToSimplifiedVerts[aTriangle.GlobalTriangle];
                TrianglesBuffer[i] = aTriangle;

                Triangles bTriangle = TrianglesBuffer[i + 1];
                bTriangle.SimplifiedTriangle = VerticesToSimplifiedVerts[bTriangle.GlobalTriangle];
                TrianglesBuffer[i + 1] = bTriangle;


                Triangles cTriangle = TrianglesBuffer[i + 2];
                cTriangle.SimplifiedTriangle = VerticesToSimplifiedVerts[cTriangle.GlobalTriangle];
                TrianglesBuffer[i + 2] = cTriangle;
                SimplifiedTrianglesBuffer.Add(new SimplifiedTriangles
                {
                    TriangleIndex = Triangles,
                    SimpVert1Index = aTriangle.SimplifiedTriangle,
                    SimpVert2Index = bTriangle.SimplifiedTriangle,
                    SimpVert3Index = cTriangle.SimplifiedTriangle
                });
                Triangles++;
            }
            #endregion

            #region Edge Data Preperation
            SimplifiedVerticesBuffer = entityManager.GetBuffer<SimplifiedVertices>(MeshManagerReal);
            DynamicBuffer<SimplifiedEdges> SimplifiedEdgesBuffer = entityManager.GetBuffer<SimplifiedEdges>(MeshManagerReal);
            NativeHashSet<Vector3> Edgepos = new NativeHashSet<Vector3>(SimplifiedVerticesBuffer.Length * 2, Allocator.Temp);
            for (int i = 0; i < TrianglesBuffer.Length; i += 3)
            {
                NativeArray<int> EdgeCollection = new NativeArray<int>(new int[] { TrianglesBuffer[i].SimplifiedTriangle, TrianglesBuffer[i + 1].SimplifiedTriangle }, Allocator.Temp);
                Vector3 TempV3 = (SimplifiedVerticesBuffer[EdgeCollection[0]].Vertex + SimplifiedVerticesBuffer[EdgeCollection[1]].Vertex) / 2;
                if (!Edgepos.Contains(TempV3))
                {
                    Edgepos.Add(TempV3);
                    SimplifiedEdgesBuffer.Add(new SimplifiedEdges { EdgeIndex = SimplifiedEdgesBuffer.Length, SimpVert1Index = EdgeCollection[0], SimpVert2Index = EdgeCollection[1], EdgePosition = TempV3 });
                }
                EdgeCollection.Dispose();
                EdgeCollection = new NativeArray<int>(new int[] { TrianglesBuffer[i + 1].SimplifiedTriangle, TrianglesBuffer[i + 2].SimplifiedTriangle }, Allocator.Temp);
                TempV3 = (SimplifiedVerticesBuffer[EdgeCollection[0]].Vertex + SimplifiedVerticesBuffer[EdgeCollection[1]].Vertex) / 2;
                if (!Edgepos.Contains(TempV3))
                {
                    Edgepos.Add(TempV3);
                    SimplifiedEdgesBuffer.Add(new SimplifiedEdges { EdgeIndex = SimplifiedEdgesBuffer.Length, SimpVert1Index = EdgeCollection[0], SimpVert2Index = EdgeCollection[1], EdgePosition = TempV3 });
                }
                EdgeCollection.Dispose();
                EdgeCollection = new NativeArray<int>(new int[] { TrianglesBuffer[i + 2].SimplifiedTriangle, TrianglesBuffer[i].SimplifiedTriangle }, Allocator.Temp);
                TempV3 = (SimplifiedVerticesBuffer[EdgeCollection[0]].Vertex + SimplifiedVerticesBuffer[EdgeCollection[1]].Vertex) / 2;
                if (!Edgepos.Contains(TempV3))
                {
                    Edgepos.Add(TempV3);
                    SimplifiedEdgesBuffer.Add(new SimplifiedEdges { EdgeIndex = SimplifiedEdgesBuffer.Length, SimpVert1Index = EdgeCollection[0], SimpVert2Index = EdgeCollection[1], EdgePosition = TempV3 });
                }

                EdgeCollection.Dispose();
            }
            #endregion

            #region Vertex Finalising
            SimplifiedVerticesBuffer = entityManager.GetBuffer<SimplifiedVertices>(MeshManagerReal);
            length = SimplifiedVerticesBuffer.Length;
            NativeArray<Entity> VertexEntities = entityManager.CreateEntity(VertDataStoreArch, length, Allocator.Temp);

            entityManager.SetSharedComponentData(GetEntityQuery(typeof(VertDataStore)), new RenderMesh() { mesh = MeshManager.staticHandleMesh, material = MeshManager.staticHandleMaterial });
            //NativeHashMap<Entity, int> simpVertIndexMap = new NativeHashMap<Entity, int>(length, Allocator.TempJob);
            for (int i = 0; i < length; i++)
            {
                SimplifiedVerticesBuffer = entityManager.GetBuffer<SimplifiedVertices>(MeshManagerReal);
                SimplifiedVertices simplifiedVertex = SimplifiedVerticesBuffer[i];
                simplifiedVertex.VertEntity = VertexEntities[i];
                SimplifiedVerticesBuffer[i] = simplifiedVertex;
                VertDataStore CurrentVDS = entityManager.GetComponentData<VertDataStore>(VertexEntities[i]);
                CurrentVDS.SimplifiedVertexIndex = i;
                entityManager.SetComponentData(VertexEntities[i], CurrentVDS);
                entityManager.SetComponentData(VertexEntities[i], new Parent() { Value = MeshManagerReal });
                DynamicBuffer<Child> children = entityManager.GetBuffer<Child>(MeshManagerReal);
                children.Add(new Child { Value = VertexEntities[i] });
//                simpVertIndexMap.Add(VertexEntities[i], i);


                SimplifiedVerticesBuffer = entityManager.GetBuffer<SimplifiedVertices>(MeshManagerReal);
                NativeArray<VDRelatedVert> vDRelatedVerts = new NativeArray<VDRelatedVert>(relatedVertsDict[SimplifiedVerticesBuffer[i].Vertex].ToArray(), Allocator.Temp);
                DynamicBuffer<VDRelatedVert> RelatedVertsBuffer = entityManager.AddBuffer<VDRelatedVert>(VertexEntities[i]);
                RelatedVertsBuffer.CopyFrom(vDRelatedVerts);
                vDRelatedVerts.Dispose();
            }
            VertexEntities.Dispose();
            NativeArray<SimplifiedVertices> SimpVertArrayFinal = entityManager.GetBuffer<SimplifiedVertices>(MeshManagerReal).ToNativeArray(Allocator.Temp);
            
            ///SpawnVertexes vertexJob = new SpawnVertexes
            ///{
            ///    VDSTypeHandle = this.GetComponentTypeHandle<VertDataStore>(),
            ///    colliderTypeHandle = this.GetComponentTypeHandle<PhysicsCollider>(),
            ///    translationTypeHandle = this.GetComponentTypeHandle<Translation>(),
            ///    vertexTypeHandle = this.GetEntityTypeHandle(),
            ///    SimplifiedVertices = SimplifiedVerticesBuffer,
            ///    entityCommandBuffer = StartUpCommandBuffer.AsParallelWriter(),
            ///    MeshManagerEntity = MeshManagerReal,
            ///    simpVertIndexMap = simpVertIndexMap
            ///};
            ///
            ///EntityQuery vertexQueryForJob = entityManager.CreateEntityQuery(VertDataStoreQuery);
            ///JobHandle VertexJobHandle = simpVertIndexMap.Dispose(vertexJob.ScheduleParallel(vertexQueryForJob, 1));
            ///VertexJobHandle.Complete();
            ///endSimulationEntityCommandBufferSystem.AddJobHandleForProducer(VertexJobHandle);
            #endregion

            #region Edge Finalising

            SimplifiedEdgesBuffer = entityManager.GetBuffer<SimplifiedEdges>(MeshManagerReal);
            length = SimplifiedEdgesBuffer.Length;
            NativeArray<Entity> EdgeEntities = entityManager.CreateEntity(EdgeDataStoreArch, length, Allocator.Temp);
            entityManager.SetSharedComponentData(GetEntityQuery(typeof(EdgeDataStore)), new RenderMesh() { mesh = MeshManager.staticHandleMesh, material = MeshManager.staticHandleMaterial });
            //NativeHashMap<Entity, int> simpEdgeIndexMap = new NativeHashMap<Entity, int>(length, Allocator.TempJob);
            for (int i = 0; i < length; i++)
            {
                SimplifiedEdgesBuffer = entityManager.GetBuffer<SimplifiedEdges>(MeshManagerReal);
                SimplifiedEdges SimplifiedEdge = SimplifiedEdgesBuffer[i];
                SimplifiedEdge.EdgeEntity = EdgeEntities[i];
                SimplifiedEdge.SimpVert1Entity = SimpVertArrayFinal[SimplifiedEdge.SimpVert1Index].VertEntity;
                SimplifiedEdge.SimpVert2Entity = SimpVertArrayFinal[SimplifiedEdge.SimpVert2Index].VertEntity;
                SimplifiedEdgesBuffer[i] = SimplifiedEdge;

                DynamicBuffer<Child> children = entityManager.GetBuffer<Child>(MeshManagerReal);
                children.Add(new Child { Value = EdgeEntities[i] });

                EdgeDataStore CurrentEDS = entityManager.GetComponentData<EdgeDataStore>(EdgeEntities[i]);
                CurrentEDS.EdgeIndex = i;
                entityManager.SetComponentData(EdgeEntities[i], CurrentEDS);
                entityManager.SetComponentData(EdgeEntities[i], new Parent() { Value = MeshManagerReal });
                //simpEdgeIndexMap.Add(EdgeEntities[i], i);

            }
            EdgeEntities.Dispose();

            ///SpawnEdges edgeJob = new SpawnEdges
            ///{
            ///    EDSTypeHandle = this.GetComponentTypeHandle<EdgeDataStore>(),
            ///    colliderTypeHandle = this.GetComponentTypeHandle<PhysicsCollider>(),
            ///    translationTypeHandle = this.GetComponentTypeHandle<Translation>(),
            ///    scaleTypeHandle = this.GetComponentTypeHandle<NonUniformScale>(),
            ///    RotationTypeHandle = this.GetComponentTypeHandle<Rotation>(),
            ///    edgeTypeHandle = this.GetEntityTypeHandle(),
            ///    SimplifiedVertices = entityManager.GetBuffer<SimplifiedVertices>(MeshManagerReal),
            ///    entityCommandBuffer = StartUpCommandBuffer.AsParallelWriter(),
            ///    MeshManagerEntity = MeshManagerReal,
            ///    SimplifiedEdges = entityManager.GetBuffer<SimplifiedEdges>(MeshManagerReal),
            ///    simpEdgeIndexMap = simpEdgeIndexMap
            ///};
            ///
            ///EntityQuery edgeQueryForJob = entityManager.CreateEntityQuery(EdgeDataStoreQuery);
            ///JobHandle EdgeJobHandle = simpEdgeIndexMap.Dispose(edgeJob.ScheduleParallel(edgeQueryForJob, 1, VertexJobHandle));
            ///endSimulationEntityCommandBufferSystem.AddJobHandleForProducer(EdgeJobHandle);
            ///EdgeJobHandle.Complete();
            #endregion

            #region Triangle Finalising
            SimplifiedTrianglesBuffer = entityManager.GetBuffer<SimplifiedTriangles>(MeshManagerReal);
            length = SimplifiedTrianglesBuffer.Length;
            NativeArray<Entity> TriangleEntities = entityManager.CreateEntity(TriangleDataStoreArch, length, Allocator.Temp);
            //NativeHashMap<Entity, int> simpTriangleIndexMap = new NativeHashMap<Entity, int>(length, Allocator.TempJob);
            for (int i = 0; i < length; i++)
            {
                SimplifiedTrianglesBuffer = entityManager.GetBuffer<SimplifiedTriangles>(MeshManagerReal);
                SimplifiedTriangles current = SimplifiedTrianglesBuffer[i];
                current.TriangleEntity = TriangleEntities[i];
                current.SimpVert1Entity = SimpVertArrayFinal[current.SimpVert1Index].VertEntity;
                current.SimpVert2Entity = SimpVertArrayFinal[current.SimpVert2Index].VertEntity;
                current.SimpVert3Entity = SimpVertArrayFinal[current.SimpVert3Index].VertEntity;
                SimplifiedTrianglesBuffer[i] = current;
                DynamicBuffer<Child> children = entityManager.GetBuffer<Child>(MeshManagerReal);
                children.Add(new Child { Value = TriangleEntities[i] });
                TriangleDataStore CurrentTDS = entityManager.GetComponentData<TriangleDataStore>(TriangleEntities[i]);
                CurrentTDS.TriangleIndex = i;
                entityManager.SetComponentData(TriangleEntities[i], CurrentTDS);
                entityManager.SetComponentData(TriangleEntities[i], new Parent() { Value = MeshManagerReal });
                //simpTriangleIndexMap.Add(TriangleEntities[i], i);
            }
            TriangleEntities.Dispose();
            SimpVertArrayFinal.Dispose();

            ///SpawnTriangles triangleJob = new SpawnTriangles
            ///{
            ///    TDSTypeHandle = this.GetComponentTypeHandle<TriangleDataStore>(),
            ///    triangleTypeHandle = this.GetEntityTypeHandle(),
            ///    SimplifiedVertices = entityManager.GetBuffer<SimplifiedVertices>(MeshManagerReal),
            ///    entityCommandBuffer = StartUpCommandBuffer.AsParallelWriter(),
            ///    MeshManagerEntity = MeshManagerReal,
            ///    SimplifiedTriangles = entityManager.GetBuffer<SimplifiedTriangles>(MeshManagerReal),
            ///    simpTriangleIndexMap = simpTriangleIndexMap
            ///};
            ///
            ///EntityQuery triangleQueryForJob = entityManager.CreateEntityQuery(TriangleDataStoreQuery);
            ///JobHandle TriangleJobHandle = simpTriangleIndexMap.Dispose(triangleJob.ScheduleParallel(triangleQueryForJob, 1, EdgeJobHandle));
            ///endSimulationEntityCommandBufferSystem.AddJobHandleForProducer(TriangleJobHandle);
            ///TriangleJobHandle.Complete();

            #endregion

            #region Jobs

            //SpawnVertexes vertexJob = new SpawnVertexes
            //{
            //    VDSTypeHandle = this.GetComponentTypeHandle<VertDataStore>(),
            //    //colliderTypeHandle = this.GetComponentTypeHandle<PhysicsCollider>(),
            //    translationTypeHandle = this.GetComponentTypeHandle<Translation>(),
            //    vertexTypeHandle = this.GetEntityTypeHandle(),
            //    SimplifiedVertices = entityManager.GetBuffer<SimplifiedVertices>(MeshManagerReal),
            //    entityCommandBuffer = StartUpCommandBuffer.AsParallelWriter(),
            //    MeshManagerEntity = MeshManagerReal,
            //    //simpVertIndexMap = simpVertIndexMap
            //};
            //
            //EntityQuery vertexQueryForJob = entityManager.CreateEntityQuery(VertDataStoreQuery);
            ////JobHandle VertexJobHandle = simpVertIndexMap.Dispose(vertexJob.ScheduleParallel(vertexQueryForJob, 1));
            //JobHandle VertexJobHandle = vertexJob.ScheduleParallel(vertexQueryForJob, 1);
            //endSimulationEntityCommandBufferSystem.AddJobHandleForProducer(VertexJobHandle);
            //VertexJobHandle.Complete();
            //
            //SpawnEdges edgeJob = new SpawnEdges
            //{
            //    EDSTypeHandle = this.GetComponentTypeHandle<EdgeDataStore>(),
            //    //colliderTypeHandle = this.GetComponentTypeHandle<PhysicsCollider>(),
            //    translationTypeHandle = this.GetComponentTypeHandle<Translation>(),
            //    scaleTypeHandle = this.GetComponentTypeHandle<NonUniformScale>(),
            //    RotationTypeHandle = this.GetComponentTypeHandle<Rotation>(),
            //    edgeTypeHandle = this.GetEntityTypeHandle(),
            //    SimplifiedVertices = entityManager.GetBuffer<SimplifiedVertices>(MeshManagerReal),
            //    entityCommandBuffer = StartUpCommandBuffer.AsParallelWriter(),
            //    MeshManagerEntity = MeshManagerReal,
            //    SimplifiedEdges = entityManager.GetBuffer<SimplifiedEdges>(MeshManagerReal),
            //    //simpEdgeIndexMap = simpEdgeIndexMap
            //};
            //
            //EntityQuery edgeQueryForJob = entityManager.CreateEntityQuery(EdgeDataStoreQuery);
            ////JobHandle EdgeJobHandle = simpEdgeIndexMap.Dispose(edgeJob.ScheduleParallel(edgeQueryForJob, 1, VertexJobHandle));
            //JobHandle EdgeJobHandle = edgeJob.ScheduleParallel(edgeQueryForJob, 1, VertexJobHandle);
            //endSimulationEntityCommandBufferSystem.AddJobHandleForProducer(EdgeJobHandle);
            //EdgeJobHandle.Complete();
            //SpawnTriangles triangleJob = new SpawnTriangles
            //{
            //    TDSTypeHandle = this.GetComponentTypeHandle<TriangleDataStore>(),
            //    triangleTypeHandle = this.GetEntityTypeHandle(),
            //    SimplifiedVertices = entityManager.GetBuffer<SimplifiedVertices>(MeshManagerReal),
            //    entityCommandBuffer = StartUpCommandBuffer.AsParallelWriter(),
            //    MeshManagerEntity = MeshManagerReal,
            //    SimplifiedTriangles = entityManager.GetBuffer<SimplifiedTriangles>(MeshManagerReal),
            //    //simpTriangleIndexMap = simpTriangleIndexMap
            //};
            //
            //EntityQuery triangleQueryForJob = entityManager.CreateEntityQuery(TriangleDataStoreQuery);
            ////JobHandle TriangleJobHandle = simpTriangleIndexMap.Dispose(triangleJob.ScheduleParallel(triangleQueryForJob, 1, EdgeJobHandle));
            //JobHandle TriangleJobHandle = triangleJob.ScheduleParallel(triangleQueryForJob, 1, EdgeJobHandle);
            //endSimulationEntityCommandBufferSystem.AddJobHandleForProducer(TriangleJobHandle);
            //TriangleJobHandle.Complete();
            #endregion

            StartUpCommandBuffer.RemoveComponent(MeshManagerReal, typeof(MeshManagerPreStart));
            StartUpCommandBuffer.AddComponent(MeshManagerReal, new ComponentTypes(typeof(MeshManagerStartVert), typeof(MeshManagerStartEdge), typeof(MeshManagerStartTri)));
            Debug.Log("MainStart() Execution Time: " + (UnityEngine.Time.realtimeSinceStartup - SetUpFirst) * 1000f + "ms");
        });
    }
}

public class StartVertexSystem : ComponentSystem
{

    private EntityManager entityManager;
    private EntityCommandBuffer StartUpCommandBuffer;
    private EndSimulationEntityCommandBufferSystem endSimulationEntityCommandBufferSystem;
    private EntityQueryDesc VertDataStoreQuery;

    protected override void OnCreate()
    {
        endSimulationEntityCommandBufferSystem = World.GetOrCreateSystem<EndSimulationEntityCommandBufferSystem>();
        entityManager = World.DefaultGameObjectInjectionWorld.EntityManager;
        VertDataStoreQuery = new EntityQueryDesc
        {
            All = new ComponentType[] { typeof(VertDataStore) }
        };
    }

    protected override void OnUpdate()
    {
        StartUpCommandBuffer = endSimulationEntityCommandBufferSystem.CreateCommandBuffer();
        Entities.WithAll<MeshManagerType, MeshManagerStartVert>().ForEach((Entity MeshManagerReal) =>
        {
            float SetUpFirst = UnityEngine.Time.realtimeSinceStartup;
            StartUpCommandBuffer.RemoveComponent(MeshManagerReal, typeof(MeshManagerStartVert));
            SpawnVertexes vertexJob = new SpawnVertexes
            {
                VDSTypeHandle = this.GetComponentTypeHandle<VertDataStore>(),
                translationTypeHandle = this.GetComponentTypeHandle<Translation>(),
                vertexTypeHandle = this.GetEntityTypeHandle(),
                SimplifiedVertices = entityManager.GetBuffer<SimplifiedVertices>(MeshManagerReal),
                entityCommandBuffer = StartUpCommandBuffer.AsParallelWriter(),
                MeshManagerEntity = MeshManagerReal,
            };
            EntityQuery vertexQueryForJob = entityManager.CreateEntityQuery(VertDataStoreQuery);
            JobHandle VertexJobHandle = vertexJob.ScheduleParallel(vertexQueryForJob, 1);
            endSimulationEntityCommandBufferSystem.AddJobHandleForProducer(VertexJobHandle);
            VertexJobHandle.Complete();
            StartUpCommandBuffer.AddComponent(MeshManagerReal, typeof(MeshManagerStartedVert));

            Debug.Log("StartVertexOnUpdate() Execution Time: " + (UnityEngine.Time.realtimeSinceStartup - SetUpFirst) * 1000f + "ms");
        });
    }
}

public class StartEdgeSystem : ComponentSystem
{

    private EntityManager entityManager;
    private EntityCommandBuffer StartUpCommandBuffer;
    private EndSimulationEntityCommandBufferSystem endSimulationEntityCommandBufferSystem;
    private EntityQueryDesc EdgeDataStoreQuery;

    protected override void OnCreate()
    {
        endSimulationEntityCommandBufferSystem = World.GetOrCreateSystem<EndSimulationEntityCommandBufferSystem>();
        entityManager = World.DefaultGameObjectInjectionWorld.EntityManager;
        EdgeDataStoreQuery = new EntityQueryDesc
        {
            All = new ComponentType[] { typeof(EdgeDataStore) }
        };
    }

    protected override void OnUpdate()
    {
        StartUpCommandBuffer = endSimulationEntityCommandBufferSystem.CreateCommandBuffer();
        Entities.WithAll<MeshManagerType, MeshManagerStartEdge>().ForEach((Entity MeshManagerReal) =>
        {
            float SetUpFirst = UnityEngine.Time.realtimeSinceStartup;
            StartUpCommandBuffer.RemoveComponent(MeshManagerReal, typeof(MeshManagerStartEdge));
            SpawnEdges edgeJob = new SpawnEdges
            {
                EDSTypeHandle = this.GetComponentTypeHandle<EdgeDataStore>(),
                translationTypeHandle = this.GetComponentTypeHandle<Translation>(),
                scaleTypeHandle = this.GetComponentTypeHandle<NonUniformScale>(),
                RotationTypeHandle = this.GetComponentTypeHandle<Rotation>(),
                edgeTypeHandle = this.GetEntityTypeHandle(),
                SimplifiedVertices = entityManager.GetBuffer<SimplifiedVertices>(MeshManagerReal),
                entityCommandBuffer = StartUpCommandBuffer.AsParallelWriter(),
                MeshManagerEntity = MeshManagerReal,
                SimplifiedEdges = entityManager.GetBuffer<SimplifiedEdges>(MeshManagerReal),
            };
            EntityQuery edgeQueryForJob = entityManager.CreateEntityQuery(EdgeDataStoreQuery);
            JobHandle EdgeJobHandle = edgeJob.ScheduleParallel(edgeQueryForJob, 1);
            endSimulationEntityCommandBufferSystem.AddJobHandleForProducer(EdgeJobHandle);
            EdgeJobHandle.Complete();
            StartUpCommandBuffer.AddComponent(MeshManagerReal, typeof(MeshManagerStartedEdge));

            Debug.Log("StartEdgeOnUpdate() Execution Time: " + (UnityEngine.Time.realtimeSinceStartup - SetUpFirst) * 1000f + "ms");
        });
    }
}

public class StartTriangleSystem : ComponentSystem
{

    private EntityManager entityManager;
    private EntityCommandBuffer StartUpCommandBuffer;
    private EndSimulationEntityCommandBufferSystem endSimulationEntityCommandBufferSystem;
    private EntityQueryDesc TriangleDataStoreQuery;

    protected override void OnCreate()
    {
        endSimulationEntityCommandBufferSystem = World.GetOrCreateSystem<EndSimulationEntityCommandBufferSystem>();
        entityManager = World.DefaultGameObjectInjectionWorld.EntityManager;
        TriangleDataStoreQuery = new EntityQueryDesc
        {
            All = new ComponentType[] { typeof(TriangleDataStore) }
        };
    }

    protected override void OnUpdate()
    {
        StartUpCommandBuffer = endSimulationEntityCommandBufferSystem.CreateCommandBuffer();
        Entities.WithAll<MeshManagerType, MeshManagerStartTri>().ForEach((Entity MeshManagerReal) =>
        {
            float SetUpFirst = UnityEngine.Time.realtimeSinceStartup;
            StartUpCommandBuffer.RemoveComponent(MeshManagerReal, typeof(MeshManagerStartTri));
            SpawnTriangles triangleJob = new SpawnTriangles
            {
                TDSTypeHandle = this.GetComponentTypeHandle<TriangleDataStore>(),
                triangleTypeHandle = this.GetEntityTypeHandle(),
                SimplifiedVertices = entityManager.GetBuffer<SimplifiedVertices>(MeshManagerReal),
                entityCommandBuffer = StartUpCommandBuffer.AsParallelWriter(),
                MeshManagerEntity = MeshManagerReal,
                SimplifiedTriangles = entityManager.GetBuffer<SimplifiedTriangles>(MeshManagerReal),
            };
            EntityQuery triangleQueryForJob = entityManager.CreateEntityQuery(TriangleDataStoreQuery);
            JobHandle TriangleJobHandle = triangleJob.ScheduleParallel(triangleQueryForJob, 1);
            endSimulationEntityCommandBufferSystem.AddJobHandleForProducer(TriangleJobHandle);
            TriangleJobHandle.Complete();
            StartUpCommandBuffer.AddComponent(MeshManagerReal, typeof(MeshManagerStartedTri));
            Debug.Log("StartTriangleOnUpdate() Execution Time: " + (UnityEngine.Time.realtimeSinceStartup - SetUpFirst) * 1000f + "ms");
        });
    }
}



[DisableAutoCreation]
public class StartUpSystem : JobComponentSystem
{
    private EntityManager entityManager;
    private EntityQueryDesc TriangleDataStoreQuery;
    private EntityQueryDesc EdgeDataStoreQuery;
    private EntityQueryDesc VertDataStoreQuery;
    protected override void OnCreate()
    {
        //this.Enabled = false;
        entityManager = World.DefaultGameObjectInjectionWorld.EntityManager;
        VertDataStoreQuery = new EntityQueryDesc
        {
            All = new ComponentType[] { typeof(VertDataStore) }
        };
        EdgeDataStoreQuery = new EntityQueryDesc
        {
            All = new ComponentType[] { typeof(EdgeDataStore) }
        };
        TriangleDataStoreQuery = new EntityQueryDesc
        {
            All = new ComponentType[] { typeof(TriangleDataStore) }
        };
    }
    protected override JobHandle OnUpdate(JobHandle inputDeps)
    {
        EndSimulationEntityCommandBufferSystem endSimulationEntityCommandBufferSystem = World.GetOrCreateSystem<EndSimulationEntityCommandBufferSystem>();
        EntityCommandBuffer StartUpCommandBuffer = endSimulationEntityCommandBufferSystem.CreateCommandBuffer();
        EntityManager localentityManager = entityManager;
        EntityQueryDesc LocalTriangleDataStoreQuery = VertDataStoreQuery;
        EntityQueryDesc LocalEdgeDataStoreQuery = EdgeDataStoreQuery;
        EntityQueryDesc LocalVertDataStoreQuery = TriangleDataStoreQuery;
        Entities.WithReadOnly(localentityManager).WithReadOnly(LocalVertDataStoreQuery).WithReadOnly(LocalEdgeDataStoreQuery).WithReadOnly(LocalTriangleDataStoreQuery).WithAll<MeshManagerType, MeshManagerStart>().ForEach((int entityInQueryIndex, Entity MeshManagerReal) =>
        {
            SpawnVertexes vertexJob = new SpawnVertexes
            {
                VDSTypeHandle = this.GetComponentTypeHandle<VertDataStore>(),
                //colliderTypeHandle = this.GetComponentTypeHandle<PhysicsCollider>(),
                translationTypeHandle = this.GetComponentTypeHandle<Translation>(),
                vertexTypeHandle = this.GetEntityTypeHandle(),
                SimplifiedVertices = localentityManager.GetBuffer<SimplifiedVertices>(MeshManagerReal),
                entityCommandBuffer = StartUpCommandBuffer.AsParallelWriter(),
                MeshManagerEntity = MeshManagerReal,
                //simpVertIndexMap = simpVertIndexMap
            };

            EntityQuery vertexQueryForJob = localentityManager.CreateEntityQuery(LocalVertDataStoreQuery);
            //JobHandle VertexJobHandle = simpVertIndexMap.Dispose(vertexJob.ScheduleParallel(vertexQueryForJob, 1));
            JobHandle VertexJobHandle = vertexJob.ScheduleParallel(vertexQueryForJob, 1, inputDeps);
            endSimulationEntityCommandBufferSystem.AddJobHandleForProducer(VertexJobHandle);

            VertexJobHandle.Complete();
            SpawnEdges edgeJob = new SpawnEdges
            {
                EDSTypeHandle = this.GetComponentTypeHandle<EdgeDataStore>(),
                //colliderTypeHandle = this.GetComponentTypeHandle<PhysicsCollider>(),
                translationTypeHandle = this.GetComponentTypeHandle<Translation>(),
                scaleTypeHandle = this.GetComponentTypeHandle<NonUniformScale>(),
                RotationTypeHandle = this.GetComponentTypeHandle<Rotation>(),
                edgeTypeHandle = this.GetEntityTypeHandle(),
                SimplifiedVertices = localentityManager.GetBuffer<SimplifiedVertices>(MeshManagerReal),
                entityCommandBuffer = StartUpCommandBuffer.AsParallelWriter(),
                MeshManagerEntity = MeshManagerReal,
                SimplifiedEdges = localentityManager.GetBuffer<SimplifiedEdges>(MeshManagerReal),
                //simpEdgeIndexMap = simpEdgeIndexMap
            };

            EntityQuery edgeQueryForJob = localentityManager.CreateEntityQuery(LocalEdgeDataStoreQuery);
            //JobHandle EdgeJobHandle = simpEdgeIndexMap.Dispose(edgeJob.ScheduleParallel(edgeQueryForJob, 1, VertexJobHandle));
            JobHandle EdgeJobHandle = edgeJob.ScheduleParallel(edgeQueryForJob, 1, VertexJobHandle);
            endSimulationEntityCommandBufferSystem.AddJobHandleForProducer(EdgeJobHandle);

            EdgeJobHandle.Complete();
            SpawnTriangles triangleJob = new SpawnTriangles
            {
                TDSTypeHandle = this.GetComponentTypeHandle<TriangleDataStore>(),
                triangleTypeHandle = this.GetEntityTypeHandle(),
                SimplifiedVertices = localentityManager.GetBuffer<SimplifiedVertices>(MeshManagerReal),
                entityCommandBuffer = StartUpCommandBuffer.AsParallelWriter(),
                MeshManagerEntity = MeshManagerReal,
                SimplifiedTriangles = localentityManager.GetBuffer<SimplifiedTriangles>(MeshManagerReal),
                //simpTriangleIndexMap = simpTriangleIndexMap
            };

            EntityQuery triangleQueryForJob = localentityManager.CreateEntityQuery(LocalTriangleDataStoreQuery);
            //JobHandle TriangleJobHandle = simpTriangleIndexMap.Dispose(triangleJob.ScheduleParallel(triangleQueryForJob, 1, EdgeJobHandle));
            JobHandle TriangleJobHandle = triangleJob.ScheduleParallel(triangleQueryForJob, 1, EdgeJobHandle);
            endSimulationEntityCommandBufferSystem.AddJobHandleForProducer(TriangleJobHandle);
            TriangleJobHandle.Complete();
            StartUpCommandBuffer.RemoveComponent(MeshManagerReal, typeof(MeshManagerStart));
            StartUpCommandBuffer.AddComponent(MeshManagerReal, typeof(MeshManagerStarted));

        }).WithoutBurst().Run();
        return inputDeps;
    }
}

[BurstCompile]
public struct SpawnTriangles : IJobEntityBatch
{
    [ReadOnly]
    public EntityTypeHandle triangleTypeHandle;
    [ReadOnly]
    public Entity MeshManagerEntity;
    [ReadOnly]
    public DynamicBuffer<SimplifiedVertices> SimplifiedVertices;
    [ReadOnly]
    public DynamicBuffer<SimplifiedTriangles> SimplifiedTriangles;
    //[ReadOnly]
    //public NativeHashMap<Entity, int> simpTriangleIndexMap;

    public ComponentTypeHandle<TriangleDataStore> TDSTypeHandle;
    public EntityCommandBuffer.ParallelWriter entityCommandBuffer;
    // A Dynamic buffer is a special type of list for DOTS
    public void Execute(ArchetypeChunk batchInChunk, int batchIndex)
    {
        // Native arrays are arrays for DOTS
        NativeArray<TriangleDataStore> TDSs = batchInChunk.GetNativeArray(TDSTypeHandle);
        NativeArray<Entity> entities = batchInChunk.GetNativeArray(triangleTypeHandle);
        for (int i = 0; i < batchInChunk.Count; i++)
        {
            int index = TDSs[i].TriangleIndex;
            NativeArray<float3> vertices = new NativeArray<float3>(3, Allocator.Temp);
            vertices[0] = SimplifiedVertices[SimplifiedTriangles[index].SimpVert1Index].Vertex;
            vertices[1] = SimplifiedVertices[SimplifiedTriangles[index].SimpVert2Index].Vertex;
            vertices[2] = SimplifiedVertices[SimplifiedTriangles[index].SimpVert3Index].Vertex;
            NativeArray<int3> triangles = new NativeArray<int3>(1, Allocator.Temp);
            triangles[0] = new int3(0, 1, 2);
            
            entityCommandBuffer.AddComponent(batchIndex, entities[i], new PhysicsCollider()
            {
                Value = Unity.Physics.MeshCollider.Create(vertices, triangles)
            });
            entityCommandBuffer.AddBuffer<TDConnectedEdges>(batchIndex, entities[i]);
        }
    }
}

[BurstCompile]
public struct SpawnEdges : IJobEntityBatch
{
    [ReadOnly]
    public EntityTypeHandle edgeTypeHandle;
    [ReadOnly]
    public Entity MeshManagerEntity;
    [ReadOnly]
    public DynamicBuffer<SimplifiedVertices> SimplifiedVertices;
    [ReadOnly]
    public DynamicBuffer<SimplifiedEdges> SimplifiedEdges;
    //[ReadOnly]
    //public NativeHashMap<Entity, int> simpEdgeIndexMap;


    public ComponentTypeHandle<EdgeDataStore> EDSTypeHandle;
    //public ComponentTypeHandle<PhysicsCollider> colliderTypeHandle;
    public ComponentTypeHandle<Translation> translationTypeHandle;
    public ComponentTypeHandle<NonUniformScale> scaleTypeHandle;
    public ComponentTypeHandle<Rotation> RotationTypeHandle;
    public EntityCommandBuffer.ParallelWriter entityCommandBuffer;
    // A Dynamic buffer is a special type of list for DOTS
    public void Execute(ArchetypeChunk batchInChunk, int batchIndex)
    {
        // Native arrays are arrays for DOTS
        NativeArray<EdgeDataStore> EDSs = batchInChunk.GetNativeArray(EDSTypeHandle);
        //NativeArray<EdgeDataStore> EDSWrites = batchInChunk.GetNativeArray(EDSTypeHandle);
        NativeArray<Translation> translations = batchInChunk.GetNativeArray(translationTypeHandle);
        //NativeArray<PhysicsCollider> colliders = batchInChunk.GetNativeArray(colliderTypeHandle);
        NativeArray<NonUniformScale> scales = batchInChunk.GetNativeArray(scaleTypeHandle);
        NativeArray<Rotation> rotations = batchInChunk.GetNativeArray(RotationTypeHandle);
        NativeArray<Entity> entities = batchInChunk.GetNativeArray(edgeTypeHandle);
        for (int i = 0; i < batchInChunk.Count; i++)
        {
            //EdgeDataStore EDS = EDSs[i];
            int index = EDSs[i].EdgeIndex;// = simpEdgeIndexMap[entities[i]];
            //EDSWrites[i] = EDS;
            //entityCommandBuffer.SetSharedComponent(batchIndex, entities[i], new RenderMesh() { mesh = MeshManager.staticHandleMesh, material = MeshManager.staticHandleMaterial });
            
            entityCommandBuffer.AddBuffer<EDConnectedTriangles>(batchIndex, entities[i]);

            float3 A = SimplifiedVertices[SimplifiedEdges[index].SimpVert1Index].Vertex;
            float3 B = SimplifiedVertices[SimplifiedEdges[index].SimpVert2Index].Vertex;
            float3 EdgeAveragePositon = (A + B) / 2;
            quaternion rotation = quaternion.LookRotationSafe(A - EdgeAveragePositon, new float3(1, 1, 1));
            float3 size = new float3(0.025f, 0.025f, math.distance(A, B));

            translations[i] = new Translation() { Value = EdgeAveragePositon };
            scales[i] = new NonUniformScale() { Value = size };
            rotations[i] = new Rotation() { Value = rotation };
            entityCommandBuffer.AddComponent(batchIndex, entities[i], new PhysicsCollider()
            {
                Value = Unity.Physics.BoxCollider.Create(new BoxGeometry
                {
                    Center = float3.zero,
                    Orientation = quaternion.identity,
                    Size = size
                })
            });
            //if (colliders[i].IsValid)
            //{
            //    colliders[i].Value.Dispose();
            //}
            //colliders[i] = new PhysicsCollider()
            //{
            //    Value = Unity.Physics.BoxCollider.Create(new BoxGeometry
            //    {
            //        Center = float3.zero,
            //        Orientation = quaternion.identity,
            //        Size = size
            //    })
            //};
        }
    }
}

[BurstCompile]
public struct SpawnVertexes : IJobEntityBatch
{
    [ReadOnly]
    public EntityTypeHandle vertexTypeHandle;
    [ReadOnly]
    public Entity MeshManagerEntity;
    [ReadOnly]
    public DynamicBuffer<SimplifiedVertices> SimplifiedVertices;
    //[ReadOnly]
    //public NativeHashMap<Entity, int> simpVertIndexMap;

    public ComponentTypeHandle<VertDataStore> VDSTypeHandle;
    //public ComponentTypeHandle<PhysicsCollider> colliderTypeHandle;
    public ComponentTypeHandle<Translation> translationTypeHandle;
    public EntityCommandBuffer.ParallelWriter entityCommandBuffer;
    // A Dynamic buffer is a special type of list for DOTS
    public void Execute(ArchetypeChunk batchInChunk, int batchIndex)
    {
        // Native arrays are arrays for DOTS
        NativeArray<VertDataStore>.ReadOnly VDSs = batchInChunk.GetNativeArray(VDSTypeHandle).AsReadOnly();        
        NativeArray<Translation> translations = batchInChunk.GetNativeArray(translationTypeHandle);
        //NativeArray<PhysicsCollider> colliders = batchInChunk.GetNativeArray(colliderTypeHandle);
        NativeArray<Entity> entities = batchInChunk.GetNativeArray(vertexTypeHandle);
        for (int i = 0; i < batchInChunk.Count; i++)
        {
            //entityCommandBuffer.SetSharedComponent(batchIndex, entities[i], new RenderMesh() { mesh = MeshManager.staticHandleMesh, material = MeshManager.staticHandleMaterial });
            
            entityCommandBuffer.SetComponent(batchIndex, entities[i], new Scale() { Value = 0.1f });
            
            entityCommandBuffer.AddBuffer<VDConnectedEdge>(batchIndex, entities[i]);
            entityCommandBuffer.AddBuffer<VDConnectedTriangle>(batchIndex, entities[i]);
            entityCommandBuffer.AddBuffer<VDConnectedVert>(batchIndex, entities[i]).Add(new VDConnectedVert { ConnectedVertIndex = VDSs[i].SimplifiedVertexIndex, ConnectedVertEntity = entities[i] });
            
            //VDS.SimplifiedVertexIndex = simpVertIndexMap[entities[i]];
            
            float3 newTranslation = SimplifiedVertices[VDSs[i].SimplifiedVertexIndex].Vertex;
            translations[i] = new Translation() { Value = newTranslation };
            entityCommandBuffer.AddComponent(batchIndex, entities[i], new PhysicsCollider()
            {
                Value = Unity.Physics.BoxCollider.Create(new BoxGeometry
                {
                    Center = float3.zero,
                    Orientation = quaternion.identity,
                    Size = new float3(0.1f)
                })
            });
            //if (colliders[i].IsValid)
            //{
            //    colliders[i].Value.Dispose();
            //}
            //colliders[i] = new PhysicsCollider()
            //{
            //    Value = Unity.Physics.BoxCollider.Create(new BoxGeometry
            //    {
            //        Center = float3.zero,
            //        Orientation = quaternion.identity,
            //        Size = new float3(0.1f)
            //    })
            //};
        }
    }
}

#endregion

#region UpdateSystems
[BurstCompile]
public class UpdateBothJobSystem : JobComponentSystem
{
    private EntityQueryDesc MeshManagerQuery;
    private EntityQueryDesc VertDataStoreQuery;
    private EntityQueryDesc EdgeDataStoreQuery;
    private EntityQueryDesc TriangleDataStoreQuery;
    private EntityQueryDesc UpdateComponentQuery;
    private EndSimulationEntityCommandBufferSystem endSimulationEntityCommandBufferSystem;
    protected override void OnCreate()
    {
        MeshManagerQuery = new EntityQueryDesc { All = new ComponentType[] { typeof(MeshManagerType) } };
        VertDataStoreQuery = new EntityQueryDesc 
        { 
            All = new ComponentType[] 
            { 
                typeof(VertDataStore), 
                typeof(Translation), 
                typeof(Scale), 
                typeof(PhysicsCollider), 
                typeof(UpdateVisualComponent) 
            },
            None = new ComponentType[]
            {
                typeof(EdgeDataStore),
                typeof(TriangleDataStore)
            }
        };
        EdgeDataStoreQuery = new EntityQueryDesc 
        { 
            All = new ComponentType[] 
            { 
                typeof(EdgeDataStore), 
                typeof(Translation),
                typeof(Rotation), 
                typeof(NonUniformScale), 
                typeof(PhysicsCollider), 
                typeof(UpdateVisualComponent) 
            },
            None = new ComponentType[]
            {
                typeof(VertDataStore),
                typeof(TriangleDataStore)
            }
        };
        TriangleDataStoreQuery = new EntityQueryDesc
        {
            All = new ComponentType[]
            {
                typeof(TriangleDataStore),
                typeof(Translation),
                typeof(Rotation),
                typeof(PhysicsCollider),
                typeof(UpdateVisualComponent)
            },
            None = new ComponentType[]
            {
                typeof(EdgeDataStore),
                typeof(VertDataStore)
            }
        };
        UpdateComponentQuery = new EntityQueryDesc { All = new ComponentType[] { typeof(UpdateVisualComponent) } };
        endSimulationEntityCommandBufferSystem = World.GetOrCreateSystem<EndSimulationEntityCommandBufferSystem>();
    }

    protected override JobHandle OnUpdate(JobHandle inputDeps)
    {
        BufferFromEntity<SimplifiedVertices> SimpVertFromEntity = GetBufferFromEntity<SimplifiedVertices>();
        BufferFromEntity<SimplifiedEdges> SimpEdgeFromEntity = GetBufferFromEntity<SimplifiedEdges>();
        BufferFromEntity<SimplifiedTriangles> SimpTriangleFromEntity = GetBufferFromEntity<SimplifiedTriangles>();
        EntityQuery MeshManager = GetEntityQuery(MeshManagerQuery);
        NativeArray<Entity> MeshManagerEntity = MeshManager.ToEntityArray(Allocator.TempJob);
        DynamicBuffer<SimplifiedVertices> simplifiedVertices = SimpVertFromEntity[MeshManagerEntity[0]];
        DynamicBuffer<SimplifiedEdges> simplifiedEdges = SimpEdgeFromEntity[MeshManagerEntity[0]];
        DynamicBuffer<SimplifiedTriangles> simplifiedTriangles = SimpTriangleFromEntity[MeshManagerEntity[0]];
        MeshManagerEntity.Dispose();

        UpdateVertexVisuals vertexJob = new UpdateVertexVisuals
        {
            VDSTypeHandle = this.GetComponentTypeHandle<VertDataStore>(true),
            translationTypeHandle = this.GetComponentTypeHandle<Translation>(),
            colliderTypeHandle = this.GetComponentTypeHandle<PhysicsCollider>(),
            scaleTypeHandle = this.GetComponentTypeHandle<Scale>(),
            SimplifiedVertices = simplifiedVertices
        };
        EntityQuery VertQuery = GetEntityQuery(VertDataStoreQuery);
        JobHandle vertexHandle = vertexJob.ScheduleParallel(VertQuery, 1, inputDeps);
        EntityQuery EdgeQuery = GetEntityQuery(EdgeDataStoreQuery);
        UpdateEdgeVisuals edgeJob = new UpdateEdgeVisuals
        {
            EDSTypeHandle = this.GetComponentTypeHandle<EdgeDataStore>(true),
            translationTypeHandle = this.GetComponentTypeHandle<Translation>(),
            scaleTypeHandle = this.GetComponentTypeHandle<NonUniformScale>(),
            colliderTypeHandle = this.GetComponentTypeHandle<PhysicsCollider>(),
            RotationTypeHandle = this.GetComponentTypeHandle<Rotation>(),
            SimplifiedVertices = simplifiedVertices,
            SimplifiedEdges = simplifiedEdges
        
        };
        JobHandle edgeHandle = edgeJob.ScheduleParallel(EdgeQuery, 1,vertexHandle);
        EntityQuery triangleQueryForJob = GetEntityQuery(TriangleDataStoreQuery);
        UpdateTriangleVisuals triangleJob = new UpdateTriangleVisuals
        {
            TDSTypeHandle = this.GetComponentTypeHandle<TriangleDataStore>(),
            colliderTypeHandle = this.GetComponentTypeHandle<PhysicsCollider>(),
            SimplifiedVertices = simplifiedVertices,
            SimplifiedTriangles = simplifiedTriangles
        };
        JobHandle triangleHandle = triangleJob.ScheduleParallel(triangleQueryForJob,1 ,edgeHandle);
        EntityQuery removeQueryForJob = GetEntityQuery(UpdateComponentQuery);
        endSimulationEntityCommandBufferSystem.CreateCommandBuffer().RemoveComponentForEntityQuery(removeQueryForJob, typeof(UpdateVisualComponent));
        return triangleHandle;
    }
}

[BurstCompile]
public struct UpdateVertexVisuals : IJobEntityBatch
{
    [ReadOnly]
    public DynamicBuffer<SimplifiedVertices> SimplifiedVertices;
    [ReadOnly]
    public ComponentTypeHandle<VertDataStore> VDSTypeHandle;
    [ReadOnly]
    public ComponentTypeHandle<Scale> scaleTypeHandle;

    public ComponentTypeHandle<PhysicsCollider> colliderTypeHandle;
    public ComponentTypeHandle<Translation> translationTypeHandle;
    // A Dynamic buffer is a special type of list for DOTS
    public void Execute(ArchetypeChunk batchInChunk, int batchIndex)
    {
        // Native arrays are arrays for DOTS
        NativeArray<VertDataStore> VDSs = batchInChunk.GetNativeArray(VDSTypeHandle);
        NativeArray<Translation> translations = batchInChunk.GetNativeArray(translationTypeHandle);
        NativeArray<PhysicsCollider> colliders = batchInChunk.GetNativeArray(colliderTypeHandle);
        
        NativeArray<Scale> scales = batchInChunk.GetNativeArray(scaleTypeHandle);
        for (int i = 0; i < batchInChunk.Count; i++)
        {
            SimplifiedVertices simplifiedVertex = SimplifiedVertices[VDSs[i].SimplifiedVertexIndex];
            float3 newTranslation = simplifiedVertex.Vertex;
            translations[i] = new Translation() { Value = newTranslation };

            if (colliders[i].IsValid)
            {
                colliders[i].Value.Dispose();
            }
            colliders[i] = new PhysicsCollider()
            {
                Value = Unity.Physics.BoxCollider.Create(new BoxGeometry
                {
                    Center = float3.zero,
                    Orientation = quaternion.identity,
                    Size = new float3(scales[i].Value)
                })
            };
        }
    }
}

[BurstCompile]
public struct UpdateEdgeVisuals : IJobEntityBatch
{
    [ReadOnly]
    public DynamicBuffer<SimplifiedVertices> SimplifiedVertices;
    [ReadOnly]
    public DynamicBuffer<SimplifiedEdges> SimplifiedEdges;
    [ReadOnly]
    public ComponentTypeHandle<EdgeDataStore> EDSTypeHandle;

    public ComponentTypeHandle<PhysicsCollider> colliderTypeHandle;
    public ComponentTypeHandle<Translation> translationTypeHandle;
    public ComponentTypeHandle<NonUniformScale> scaleTypeHandle;
    public ComponentTypeHandle<Rotation> RotationTypeHandle;

    public void Execute(ArchetypeChunk batchInChunk, int batchIndex)
    {
        NativeArray<EdgeDataStore> EDSs = batchInChunk.GetNativeArray(EDSTypeHandle);
        NativeArray<Translation> translations = batchInChunk.GetNativeArray(translationTypeHandle);
        NativeArray<PhysicsCollider> colliders = batchInChunk.GetNativeArray(colliderTypeHandle);
        NativeArray<NonUniformScale> scales = batchInChunk.GetNativeArray(scaleTypeHandle);
        NativeArray<Rotation> rotations = batchInChunk.GetNativeArray(RotationTypeHandle);
        for (int i = 0; i < batchInChunk.Count; i++)
        {
            SimplifiedEdges simplifiedEdge = SimplifiedEdges[EDSs[i].EdgeIndex];
            float3 A = SimplifiedVertices[simplifiedEdge.SimpVert1Index].Vertex;
            float3 B = SimplifiedVertices[simplifiedEdge.SimpVert2Index].Vertex;
            float3 EdgeAveragePositon = (A + B) / 2;
            quaternion rotation = quaternion.LookRotationSafe(A - EdgeAveragePositon, new float3(1, 1, 1));
            float3 size = new float3(0.025f, 0.025f, math.distance(A, B));
            translations[i] = new Translation() { Value = EdgeAveragePositon };
            scales[i] = new NonUniformScale() { Value = size };
            rotations[i] = new Rotation() { Value = rotation };

            if (colliders[i].IsValid)
            {
                colliders[i].Value.Dispose();
            }
            colliders[i] = new PhysicsCollider()
            {
                Value = Unity.Physics.BoxCollider.Create(new BoxGeometry
                {
                    Center = float3.zero,
                    Orientation = quaternion.identity,
                    Size = size
                })
            };
        }
    }
}

[BurstCompile]
public struct UpdateTriangleVisuals : IJobEntityBatch
{
    [ReadOnly]
    public DynamicBuffer<SimplifiedVertices> SimplifiedVertices;    
    [ReadOnly]
    public DynamicBuffer<SimplifiedTriangles> SimplifiedTriangles;

    public ComponentTypeHandle<TriangleDataStore> TDSTypeHandle;
    public ComponentTypeHandle<PhysicsCollider> colliderTypeHandle;

    public void Execute(ArchetypeChunk batchInChunk, int batchIndex)
    {
        NativeArray<TriangleDataStore> TDSs = batchInChunk.GetNativeArray(TDSTypeHandle);
        NativeArray<PhysicsCollider> colliders = batchInChunk.GetNativeArray(colliderTypeHandle);
        for (int i = 0; i < batchInChunk.Count; i++)
        {
            TriangleDataStore TDS = TDSs[i];
            SimplifiedTriangles simplifiedTriangle = SimplifiedTriangles[TDS.TriangleIndex];
            float3 A = SimplifiedVertices[simplifiedTriangle.SimpVert1Index].Vertex;
            float3 B = SimplifiedVertices[simplifiedTriangle.SimpVert2Index].Vertex;
            float3 C = SimplifiedVertices[simplifiedTriangle.SimpVert3Index].Vertex;
            TDS.TrianglePosition = (A + B + C) / 3;
            TDSs[i] = TDS;
            NativeArray<float3> vertices = new NativeArray<float3>(3, Allocator.Temp);
            vertices[0] = A;
            vertices[1] = B;
            vertices[2] = C;
            NativeArray<int3> triangles = new NativeArray<int3>(1, Allocator.Temp);
            triangles[0] = new int3(0, 1, 2);
            if (colliders[i].IsValid)
            {
                colliders[i].Value.Dispose();
            }
            colliders[i] = new PhysicsCollider()
            {
                Value = Unity.Physics.MeshCollider.Create(vertices, triangles)
            };
        }
    }
}
#endregion