//Copyright(c) <2022> <JackyGun:twitter@konchannyan>
//Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files(the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and / or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions :
//
//The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Shader "JackyGun/LeadParticle/Ohhagari"
{
	Properties
	{
		_IOR("IOR (屈折率)", Float) = 1.5
		_Color("Color", Color) = (0.4, 0.6, 0.7, 1.0)
		_Opacity("Opacity (不透明度)", Range(0.0, 1.0)) = 0.0
		_Alpha("Alpha", Range(0.0, 1.0)) = 1.0

		_RateFresnel("RateFresnel", Float) = 0.076
		_FixedR("FixedR (反射方向取り入れ)", Range(0.0, 1.0)) = 0
		_FixedO("FixedO (屈折方向取り入れ)", Range(0.0, 1.0)) = 0
	}
	
	SubShader
	{
		Tags
		{
			"RenderType" = "Opaque"
			"DisableBatching" = "True"
			"Queue" = "Transparent+100"
		}

		Pass
		{
			Cull Off
			Blend One OneMinusSrcAlpha
		
			CGPROGRAM

#define GP_INSTANCE_SIZE 64

#pragma target 5.0

#pragma vertex mainVS
#pragma geometry mainGS
#pragma fragment mainFS

#pragma multi_compile_instancing
#pragma instancing_options procedural:vertInstancingSetup

#include "UnityCG.cginc"
#include "UnityStandardParticleInstancing.cginc"

			struct VS_IN
			{
				uint  instanceID : SV_InstanceID;
			};

			struct VS_OUT
			{
				uint  instanceID : SV_InstanceID;
			};

			struct GS_OUT
			{
				float4 vertex : SV_POSITION;
				float4 world  : WORLD;
				float3 model  : MODEL;
				float3 bone0  : BONE0;
				float3 bone1  : BONE1;
				float3 bone2  : BONE2;
				float3 bone3  : BONE3;
				float3 rotation : ROTATION;
			};

			struct FS_OUT
			{
				float4 color : SV_Target;
				float depth : SV_Depth;
			};

			VS_OUT mainVS(VS_IN In)
			{
				VS_OUT Out;
				Out.instanceID = In.instanceID;
				return Out;
			}

			[maxvertexcount(24)] [instance(1)]
			void mainGS(point VS_OUT input[1], inout TriangleStream<GS_OUT> outStream, uint gsid : SV_GSInstanceID, uint svpid : SV_PrimitiveID)
			{
				GS_OUT o = (GS_OUT)0;
				VS_OUT dm = (VS_OUT)0;

				int i, j;

				float4 bone = float4(0, 0, 10.5, 1);

				const int BONE_COUNT = 4;
				const int P_LENGTH = 32;

				const float BONE_LENGTH[BONE_COUNT] = { 0.0f, 0.25f, 0.20f, 0.05f};

				uint iid = input[0].instanceID;
				if (iid > 0) return;

				dm.instanceID = 0;
				UNITY_SETUP_INSTANCE_ID(dm);
				float4 lp4 = mul(unity_ObjectToWorld, bone);
				float3 lp = lp4.xyz / lp4.w;

				float3 pp[BONE_COUNT] =
				{
					float3(0,0,0), float3(0,0,0), float3(0,0,0), float3(0,0,0),
				};
				float3 vp[BONE_COUNT] =
				{
					float3(0,0,0), float3(0,0,0), float3(0,0,0), float3(0,0,0),
				};

				dm.instanceID = P_LENGTH - 1;
				UNITY_SETUP_INSTANCE_ID(dm);
				float4 bp4 = mul(unity_ObjectToWorld, bone);
				float3 bp = bp4.xyz / bp4.w;
				for (i = 1; i < P_LENGTH; i++)
				{
					dm.instanceID = P_LENGTH - 1 - i;
					UNITY_SETUP_INSTANCE_ID(dm);
					float4 tp4 = mul(unity_ObjectToWorld, bone);
					float3 tp = tp4.xyz / tp4.w;
					float3 tv = tp - bp;

					pp[0] = tv;
					for (j = 1; j < BONE_COUNT; j++)
					{
						float3 ttd = pp[j] - pp[j - 1] - tv;
						ttd += float3(0, -0.15f, 0);
						float3 ttv = length(ttd) > 0.0001f ? normalize(ttd) * min(BONE_LENGTH[j], length(ttd)) : 0;
						ttv = lerp(ttv, vp[j], (float)(j) / (BONE_COUNT - 1) * 0.5f);
						pp[j] = pp[j - 1] + ttv;
						vp[j] = ttv;
					}
					bp = tp;
				}

				o.model = unity_ObjectToWorld._14_24_34 / unity_ObjectToWorld._44;
				o.bone0 = lp + pp[0];
				o.bone1 = lp + pp[1];
				o.bone2 = lp + pp[2];
				o.bone3 = lp + pp[3];

				float3 v12 = normalize(o.bone1 - o.bone2);

				o.rotation = float3(atan2(v12.x, v12.y), 0, atan2(v12.y, v12.z));

				float3 pos000 = lp + float3(-1, -1, -1);
				float3 pos100 = lp + float3(+1, -1, -1);
				float3 pos010 = lp + float3(-1, +1, -1);
				float3 pos001 = lp + float3(-1, -1, +1);
				float3 pos011 = lp + float3(-1, +1, +1);
				float3 pos101 = lp + float3(+1, -1, +1);
				float3 pos110 = lp + float3(+1, +1, -1);
				float3 pos111 = lp + float3(+1, +1, +1);

				o.world = float4(pos010, 1); o.vertex = mul(UNITY_MATRIX_VP, o.world); outStream.Append(o);
				o.world = float4(pos110, 1); o.vertex = mul(UNITY_MATRIX_VP, o.world); outStream.Append(o);
				o.world = float4(pos000, 1); o.vertex = mul(UNITY_MATRIX_VP, o.world); outStream.Append(o);
				o.world = float4(pos100, 1); o.vertex = mul(UNITY_MATRIX_VP, o.world); outStream.Append(o);
				outStream.RestartStrip();

				o.world = float4(pos110, 1); o.vertex = mul(UNITY_MATRIX_VP, o.world); outStream.Append(o);
				o.world = float4(pos111, 1); o.vertex = mul(UNITY_MATRIX_VP, o.world); outStream.Append(o);
				o.world = float4(pos100, 1); o.vertex = mul(UNITY_MATRIX_VP, o.world); outStream.Append(o);
				o.world = float4(pos101, 1); o.vertex = mul(UNITY_MATRIX_VP, o.world); outStream.Append(o);
				outStream.RestartStrip();

				o.world = float4(pos111, 1); o.vertex = mul(UNITY_MATRIX_VP, o.world); outStream.Append(o);
				o.world = float4(pos011, 1); o.vertex = mul(UNITY_MATRIX_VP, o.world); outStream.Append(o);
				o.world = float4(pos101, 1); o.vertex = mul(UNITY_MATRIX_VP, o.world); outStream.Append(o);
				o.world = float4(pos001, 1); o.vertex = mul(UNITY_MATRIX_VP, o.world); outStream.Append(o);
				outStream.RestartStrip();

				o.world = float4(pos011, 1); o.vertex = mul(UNITY_MATRIX_VP, o.world); outStream.Append(o);
				o.world = float4(pos010, 1); o.vertex = mul(UNITY_MATRIX_VP, o.world); outStream.Append(o);
				o.world = float4(pos001, 1); o.vertex = mul(UNITY_MATRIX_VP, o.world); outStream.Append(o);
				o.world = float4(pos000, 1); o.vertex = mul(UNITY_MATRIX_VP, o.world); outStream.Append(o);
				outStream.RestartStrip();

				o.world = float4(pos011, 1); o.vertex = mul(UNITY_MATRIX_VP, o.world); outStream.Append(o);
				o.world = float4(pos111, 1); o.vertex = mul(UNITY_MATRIX_VP, o.world); outStream.Append(o);
				o.world = float4(pos010, 1); o.vertex = mul(UNITY_MATRIX_VP, o.world); outStream.Append(o);
				o.world = float4(pos110, 1); o.vertex = mul(UNITY_MATRIX_VP, o.world); outStream.Append(o);
				outStream.RestartStrip();

				o.world = float4(pos000, 1); o.vertex = mul(UNITY_MATRIX_VP, o.world); outStream.Append(o);
				o.world = float4(pos100, 1); o.vertex = mul(UNITY_MATRIX_VP, o.world); outStream.Append(o);
				o.world = float4(pos001, 1); o.vertex = mul(UNITY_MATRIX_VP, o.world); outStream.Append(o);
				o.world = float4(pos101, 1); o.vertex = mul(UNITY_MATRIX_VP, o.world); outStream.Append(o);
				outStream.RestartStrip();

			}

			float _IOR;
			float4 _Color;
			float _Opacity;
			float _Alpha;

			float _RateFresnel;
			float _FixedR;
			float _FixedO;

			float3x3 MatRollX(float a)
			{
				return float3x3(1, 0, 0, 0, cos(a), -sin(a), 0, sin(a), cos(a));
			}

			float3x3 MatRollY(float a)
			{
				return float3x3(cos(a), 0, sin(a), 0, 1, 0, -sin(a), 0, cos(a));
			}

			float3x3 MatRollZ(float a)
			{
				return float3x3(cos(a), -sin(a), 0, sin(a), cos(a), 0, 0, 0, 1);
			}

			float4 SampleP1(float3 i)
			{
				i = mul((float3x3)unity_WorldToObject, i);

				float4 refColor = UNITY_SAMPLE_TEXCUBE_LOD(unity_SpecCube0, i, 0);
				refColor = float4(DecodeHDR(refColor, unity_SpecCube0_HDR).rgb, 1);
				return refColor;

			}

			float FresnelApprox(float3 i, float3 n, float IOR)
			{
				float r = ((1 - (1 / IOR)) * (1 - (1 / IOR))) / ((1 + (1 / IOR)) * (1 + (1 / IOR)));
				return r + (1 - r) * pow(1 - dot(i, n), 5.0);
			}

			bool HitSphere(float3 s_pos, float s_rad, float3 v_pos, float3 v_vec, bool first, out float3 hit_pos, out float hit_dis, out float3 hit_normal)
			{
				float3 d = normalize(v_vec);
				float a = 1;
				float3 s = v_pos - s_pos;
				float b = 2 * dot(s, d);
				float c = length(s) * length(s) - s_rad * s_rad;
				float D = b * b - 4 * a * c;
				if (D < 0) {
					return false;
				}
				float t0 = (-b - sqrt(D)) / (2 * a);
				float t1 = (-b + sqrt(D)) / (2 * a);
				float tmin = min(t0, t1);
				float tmax = max(t0, t1);
				float t = first ? tmin : tmax;
				if (first && tmin < 0) {
					return false;
				}
				hit_pos = v_pos + t * d;
				hit_normal = normalize(hit_pos - s_pos);
				hit_dis = distance(v_pos, hit_pos);
				return true;
			}

			bool HitCage(float3x3 mat1, float3x3 mat2, float3 m_pos, float3 v_pos, float3 v_vec, bool first, out float3 hit_pos, out float hit_dis, out float3 hit_normal)
			{
				int i;
				float e = 0.00001f;

				float c_v_pos = v_pos;

				v_pos -= m_pos;
				v_pos = mul(mat1, v_pos);
				v_vec = mul(mat1, v_vec);

				float scale_o = 0.20f;
				float scale_i = 0.15f;

				float3 cage_positions[6] =
				{
					float3(-1, 0, 0), float3(+1, 0, 0), float3(0, -1, 0), float3(0, +1, 0), float3(0, 0, -1), float3(0, 0, +1),
				};

				bool is_hit = false;
				float b_distance = 0;

				for (i = 0; i < 6; i++)
				{
					float3 cp = cage_positions[i] * scale_o;
					float3 cn = cage_positions[i];
					float3 lvp = v_pos - cp;
					float3 dvn = dot(cn, v_vec);
					if (any(abs(dvn) < 0.000001)) { continue; }
					float tvn = -dot(lvp, cn) / dvn;
					float3 tp = lvp + v_vec * tvn;
					if (any(abs(tp) > scale_o) || all(abs(tp) < scale_i)) { continue; }

					float3 tmp_hit = tp + cp;
					float tmp_dis = distance(v_pos, tmp_hit);
					float3 tmp_normal = cn;
					
					if (!is_hit || b_distance > tmp_dis)
					{
						hit_pos = tmp_hit;
						hit_dis = tmp_dis;
						hit_normal = tmp_normal;

						is_hit = true;
						b_distance = tmp_dis;
					}
				}

				for (i = 0; i < 6; i++)
				{
					float3 cp = cage_positions[i] * scale_i;
					float3 cn = cage_positions[i];
					float3 lvp = v_pos - cp;
					float3 dvn = dot(cn, v_vec);
					if (any(abs(dvn) < 0.000001)) { continue; }
					float tvn = -dot(lvp, cn) / dvn;
					float3 tp = lvp + v_vec * tvn;
					if (any(abs(tp) > scale_o) || all(abs(tp) < scale_i)) { continue; }

					float3 tmp_hit = tp + cp;
					float tmp_dis = distance(v_pos, tmp_hit);
					float3 tmp_normal = cn;

					if (!is_hit || b_distance > tmp_dis)
					{
						hit_pos = tmp_hit;
						hit_dis = tmp_dis;
						hit_normal = -tmp_normal;

						is_hit = true;
						b_distance = tmp_dis;
					}
				}

				hit_pos = mul(mat2, hit_pos);
				hit_pos += m_pos;
				hit_normal = mul(mat2, hit_normal);
				hit_dis = distance(hit_pos, c_v_pos);
				return is_hit;
			}

			float rand(float2 p, float seed)
			{
				return frac(sin(dot(p, fixed2(12.9898, 78.233)) + seed) * 43758.5453);
			}

			float random(float2 st, float seed)
			{
				float2 p = floor(st);
				float2 f = frac(st);
				float2 u = f * f*(3.0 - 2.0*f);

				float v00 = rand(p + float2(0, 0), seed);
				float v10 = rand(p + float2(1, 0), seed);
				float v01 = rand(p + float2(0, 1), seed);
				float v11 = rand(p + float2(1, 1), seed);

				return lerp(lerp(dot(v00, f - float2(0, 0)), dot(v10, f - float2(1, 0)), u.x),
					lerp(dot(v01, f - float2(0, 1)), dot(v11, f - float2(1, 1)), u.x),
					u.y) + 0.5f;
			}

			float4 HsvColor(float3 hsv)
			{
				float4 rgb = 1;

				if (hsv.y == 0) {
					rgb.r = rgb.g = rgb.b = hsv.z;
				}
				else {
					hsv.x *= 6.0;
					float i = floor(hsv.x);
					float f = hsv.x - i;
					float aa = hsv.z * (1 - hsv.y);
					float bb = hsv.z * (1 - (hsv.y * f));
					float cc = hsv.z * (1 - (hsv.y * (1 - f)));
					if (i < 1) {
						rgb.r = hsv.z;
						rgb.g = cc;
						rgb.b = aa;
					}
					else if (i < 2) {
						rgb.r = bb;
						rgb.g = hsv.z;
						rgb.b = aa;
					}
					else if (i < 3) {
						rgb.r = aa;
						rgb.g = hsv.z;
						rgb.b = cc;
					}
					else if (i < 4) {
						rgb.r = aa;
						rgb.g = bb;
						rgb.b = hsv.z;
					}
					else if (i < 5) {
						rgb.r = cc;
						rgb.g = aa;
						rgb.b = hsv.z;
					}
					else {
						rgb.r = hsv.z;
						rgb.g = aa;
						rgb.b = bb;
					}
				}
				return rgb;
			}

			float4 PastelColor(float3 seed)
			{
				float s = random(float2(seed.x - seed.y, seed.z + seed.y) * 3 , 0);

				return HsvColor(float3(s, 0.3, 0.9));
			}

			FS_OUT mainFS(GS_OUT In, float face : VFACE)
			{
				float radius = length(unity_ObjectToWorld[0].xyz) * 0.15f;

				FS_OUT Out;
				Out.color = 0;
				Out.depth = 1;

				float3 EyePos = _WorldSpaceCameraPos;
				float3 v = normalize(In.world.xyz - EyePos);
				float3 p = face > 0 ? In.world.xyz : EyePos;
				float3 spos = In.bone3.xyz;
				float3 hit_pos = 0;
				float hit_dis = 0;
				float3 hit_normal = 0;
				float3 ipos[3];
				float3 ivec[3];
				float3 opos[3];
				float3 ovec[3];
				float3 ni[3];
				float3 no[3];
				float3 ri[3];
				float3 ro[3];

				float ior[3];
				ior[0] = _IOR * 0.99f;
				ior[1] = _IOR * 1.00f;
				ior[2] = _IOR * 1.01f;

				float3x3 mat1 = mul(MatRollX(+In.rotation.z), MatRollZ(+In.rotation.x));
				float3x3 mat2 = mul(MatRollZ(-In.rotation.x), MatRollX(-In.rotation.z));

				float3 final_hit_pos = 0;
				float3 final_hit_normal = 0;
				float final_distance = 100000000;
				for (int i = 0; i < 10; i++)
				{
					float radiusT = length(unity_ObjectToWorld[0].xyz) * 0.015f;
					if (HitSphere(lerp(In.model, In.bone0, float3(i / 10.0f, sqrt(i / 10.0f), i / 10.0f)), radiusT, p, v, true, hit_pos, hit_dis, hit_normal))
					{
						float td = distance(hit_pos, p);
						if (final_distance > td)
						{
							final_hit_pos = hit_pos;
							final_distance = td;
							Out.color = i % 2 ? float4(0.95, 0.55, 0.7, 1) : i == 0 ? float4(1, 0.1, 0.1, 1) :  float4(0.45, 0.6, 0.95, 1);
						}
					}
				}

				for (int i = 0; i < 6; i++)
				{
					float radiusT = length(unity_ObjectToWorld[0].xyz) * 0.015f;
					if (HitSphere(lerp(In.bone0, In.bone1, i / 6.0f), radiusT, p, v, true, hit_pos, hit_dis, hit_normal))
					{
						float td = distance(hit_pos, p);
						if (final_distance > td)
						{
							final_hit_pos = hit_pos;
							final_distance = td;
							Out.color = i % 2 ? float4(0.95, 0.55, 0.7, 1) : float4(0.45, 0.6, 0.95, 1);
						}
					}
				}

				if (HitCage(mat1, mat2, In.bone2, p, v, true, hit_pos, hit_dis, hit_normal))
				{
					float td = distance(hit_pos, p);
					if (final_distance > td)
					{
						final_hit_pos = hit_pos;
						final_distance = td;
						final_hit_normal = hit_normal;

						Out.color = PastelColor(final_hit_pos);
					}
				}

				if (HitSphere(spos, radius, p, v, true, hit_pos, hit_dis, hit_normal)) {

					ni[0] = hit_normal;
					ni[1] = hit_normal;
					ni[2] = hit_normal;

					ipos[0] = hit_pos;
					ipos[1] = hit_pos;
					ipos[2] = hit_pos;

					ivec[0] = normalize(refract(v, ni[0], 1 / ior[0]));
					ivec[1] = normalize(refract(v, ni[1], 1 / ior[1]));
					ivec[2] = normalize(refract(v, ni[2], 1 / ior[2]));

					HitSphere(spos, radius, ipos[0], ivec[0], false, hit_pos, hit_dis, hit_normal);

					no[0] = hit_normal;
					opos[0] = hit_pos;

					HitSphere(spos, radius, ipos[1], ivec[1], false, hit_pos, hit_dis, hit_normal);

					no[1] = hit_normal;
					opos[1] = hit_pos;

					HitSphere(spos, radius, ipos[2], ivec[2], false, hit_pos, hit_dis, hit_normal);

					no[2] = hit_normal;
					opos[2] = hit_pos;

					ro[0] = refract(ivec[0], -no[0], ior[0]);
					ro[1] = refract(ivec[1], -no[1], ior[1]);
					ro[2] = refract(ivec[2], -no[2], ior[2]);

					ri[0] = reflect(v, ni[0]);
					ri[1] = reflect(v, ni[1]);
					ri[2] = reflect(v, ni[2]);

					float4 iCol[3];
					iCol[0] = SampleP1(ri[0]);
					iCol[1] = SampleP1(ri[1]);
					iCol[2] = SampleP1(ri[2]);
					float4 oCol[3];
					oCol[0] = SampleP1(opos[0] - spos + ro[0]);
					oCol[1] = SampleP1(opos[1] - spos + ro[1]);
					oCol[2] = SampleP1(opos[2] - spos + ro[2]);

					if (HitCage(mat1, mat2, In.bone2, ipos[0], ri[0], true, hit_pos, hit_dis, hit_normal))
					{
						iCol[0] = PastelColor(hit_pos);
					}

					if (HitCage(mat1, mat2, In.bone2, ipos[1], ri[1], true, hit_pos, hit_dis, hit_normal))
					{
						iCol[1] = PastelColor(hit_pos);
					}

					if (HitCage(mat1, mat2, In.bone2, ipos[2], ri[2], true, hit_pos, hit_dis, hit_normal))
					{
						iCol[2] = PastelColor(hit_pos);
					}

					if (HitCage(mat1, mat2, In.bone2, opos[0], ro[0], true, hit_pos, hit_dis, hit_normal))
					{
						oCol[0] = PastelColor(hit_pos);
					}

					if (HitCage(mat1, mat2, In.bone2, opos[1], ro[1], true, hit_pos, hit_dis, hit_normal))
					{
						oCol[1] = PastelColor(hit_pos);
					}

					if (HitCage(mat1, mat2, In.bone2, opos[2], ro[2], true, hit_pos, hit_dis, hit_normal))
					{
						oCol[2] = PastelColor(hit_pos);
					}

					float fresnel[3];
					fresnel[0] = FresnelApprox(normalize(v), normalize(ni[0]), ior[0]);
					fresnel[1] = FresnelApprox(normalize(v), normalize(ni[1]), ior[1]);
					fresnel[2] = FresnelApprox(normalize(v), normalize(ni[2]), ior[2]);
					fresnel[0] *= _RateFresnel;
					fresnel[1] *= _RateFresnel;
					fresnel[2] *= _RateFresnel;

					fresnel[0] = saturate(fresnel[0]);
					fresnel[1] = saturate(fresnel[1]);
					fresnel[2] = saturate(fresnel[2]);

					float4 iColFinal = float4(iCol[0].r, iCol[1].g, iCol[2].b, 1);
					float4 oColFinal = float4(oCol[0].r, oCol[1].g, oCol[2].b, 1);
					oColFinal = lerp(oColFinal, _Color, _Opacity);
					float4 fCol = float4(
						lerp(iColFinal.r, oColFinal.r, fresnel[0]),
						lerp(iColFinal.g, oColFinal.g, fresnel[1]),
						lerp(iColFinal.b, oColFinal.b, fresnel[2]),
						1 
					);

					float4 Col = fCol;
					Col = lerp(Col, iColFinal, _FixedR);
					Col = lerp(Col, oColFinal, _FixedO);

					float td = distance(ipos[0], p);
					if (final_distance > td)
					{
						final_hit_pos = ipos[0];
						final_distance = td;
						Out.color = float4(Col.rgb, _Alpha);
					}
				}

				if (final_distance > 10000)
				{
					discard;
				}

				float4 sv_depth = mul(UNITY_MATRIX_VP, float4(final_hit_pos, 1));
				float h_depth = sv_depth.z / sv_depth.w;

				Out.depth = h_depth;
				return Out;
			}

			ENDCG
		}
	}
}
