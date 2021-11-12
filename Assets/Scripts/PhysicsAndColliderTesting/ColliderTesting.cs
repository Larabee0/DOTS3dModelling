using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Entities;
using Unity.Transforms;
using Unity.Collections;
using Unity.Rendering;
using Unity.Mathematics;
using System.Linq;
using Unity.Physics;

public class ColliderTesting : MonoBehaviour
{
    private EntityManager entityManager;
    private EntityArchetype MeshManagerArch;
    private Entity SpawnedMeshManagerArch;
    [SerializeField] private UnityEngine.Material[] materials;
    [SerializeField] private Mesh handleMesh;
    [SerializeField] private Mesh originalMesh;
    // Start is called before the first frame update
    void Start()
    {
        entityManager = World.DefaultGameObjectInjectionWorld.EntityManager;
        CreateEntityArchetype();
        SpawnAndSetEntity();
    }

    // Update is called once per frame
    void Update()
    {
        entityManager.SetComponentData(SpawnedMeshManagerArch, new Scale
        {
            Value = 2
        });
        entityManager.SetComponentData(SpawnedMeshManagerArch, new PhysicsCollider
        {
            Value = Unity.Physics.BoxCollider.Create(new BoxGeometry
            {
                Center = float3.zero,
                Orientation = quaternion.identity,
                Size = new float3(2, 2, 2)
            })
        });
    }
    private void SpawnAndSetEntity()
    {
        SpawnedMeshManagerArch = entityManager.CreateEntity(MeshManagerArch);

        entityManager.SetComponentData(SpawnedMeshManagerArch, new Translation
        {
            Value = new float3(transform.position)
        });
        entityManager.SetComponentData(SpawnedMeshManagerArch, new Scale
        {
            Value = 1
        });
        entityManager.SetSharedComponentData(SpawnedMeshManagerArch, new RenderMesh
        {
            mesh = originalMesh,
            material = materials[0]
        });
        BlobAssetReference<Unity.Physics.Collider> collider = Unity.Physics.BoxCollider.Create(new BoxGeometry
        {
            Center = float3.zero,
            Orientation = quaternion.identity,
            Size = new float3(1, 1, 1)
        });
        entityManager.SetComponentData(SpawnedMeshManagerArch, new PhysicsCollider { Value = collider });

    }
    private void CreateEntityArchetype()
    {
        MeshManagerArch = entityManager.CreateArchetype(
            typeof(Translation),
            typeof(Scale),
            typeof(RenderMesh),
            typeof(LocalToWorld),
            typeof(RenderBounds),
            typeof(PhysicsCollider),
            typeof(Child)
        );
    }
    private void SetScale()
    {
        
    }
    private void OnGUI()
    {
        if (GUI.Button(new Rect(10, 0, 120, 30), "Set Scale"))
        {
            SetScale();
        }
    }
}
