#ifndef SHA256_H
#define SHA256_H
/** @file sha256sum�̃v���g�^�C�v���錾����Ă��܂�
	
	@date 2005-11-20
	@auther Naoki_I
 */

typedef unsigned int u_32;

//64bit�������֌W

//long long �� LONG LONG���g���Ȃ����ł�
//NON_SUPPORT_64�̃R�����g�A�E�g���O���Ă�������
//#define NON_SUPPORT_64

#ifndef NON_SUPPORT_64
	
	#ifdef WIN_32
	typedef unsigned long long u_64;
	#else
	typedef unsigned long long u_64;
	#endif
	#define U64_HIGH(a) ((u_32)(m_MesLen>>32))
	#define U64_LOW(a) ((u_32)a&0xFFFFFFFF)

#else

/**
	�����Z�݂̂ɑΉ�����64bit����
	+, =, +=, <<(�X�g���[���o��) ���Z�q���I�[�o�[���[�h�ς�
 */
#define U64_HIGH(a) ((u_32)a.High)
#define U64_LOW(a) ((u_32)a.Low)
struct u_64
{
	u_32 High, Low;	
	u_64()
	{
		High=0;
		Low=0;
	}
	u_64(u_32 h, u_32 l)
	{
		High=h;
		Low=l;
	}
	const u_64& operator =(const u_32&n)
	{
		High=0;
		Low=n;
		return *this;
	}
	const u_64 operator +(const u_32&n)
	{
		u_64 t(High, Low);
		t+=n;
		return t;
	}
	const u_64& operator +=(const u_32&n)
	{
		if(Low>0xFFFFFFFF-n)
			High++;
		Low+=n;
		return *this;
	}
};
#endif //#ifndef NON_SUPPORT_64

/**
	SHA256Hash�����߂�N���X
	
	Push���\�b�h�ŏ������郁�b�Z�[�W���w�肷��D
	�������b�Z�[�W�̃n�b�V�������߂�ꍇ���́C
	������ɕ�����Push���\�b�h��ǂ�ł��悢�D

	�Ō��Final���\�b�h���ĂԂƃn�b�V���l��Ԃ��܂��D
	@date 2005-11-21
	@auther Naoki_I
 */
class SHA256
{
	u_32 H[8];

	u_32 W[64];
	u_32 m_WIndex;///���݌v�Z����W�̃C���f�b�N�X(�������C�P�ʂ̓o�C�g)
	u_64 m_MesLen;

	/**
		intermdediate hash value�����߂�D
		���ʂ�H[?]�ɑ�������D
	 */
	void CalcIntermediateHash();
public:
	SHA256();

	/**
		���b�Z�[�W��ǉ�����

		@param pMes �ǉ����郁�b�Z�[�W�̐擪�̃|�C���^
		@param Len	�ǉ����郁�b�Z�[�W�̒���
	 */
	void Push(unsigned char*pMes, u_32 Len);

	/**
		Push���ꂽ���b�Z�[�W����n�b�V���l�����߂�

		@param pHash ���ʂ�������z��ւ̃|�C���^(256bit�K�v�ł�)
	 */
	void Final(unsigned char* pHash);
};

#endif
