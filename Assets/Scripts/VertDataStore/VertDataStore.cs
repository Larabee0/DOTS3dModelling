using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Entities;
using Unity.Transforms;
using Unity.Mathematics;
using Unity.Collections;

public struct VertDataStore : IComponentData
{
    public int SimplifiedVertexIndex;
}
